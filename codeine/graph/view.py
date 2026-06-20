import random

from dataclasses import dataclass
from typing import Dict, Generator, List, Optional, Sequence, Tuple, Union

from codeine.graph.graph import CodonGraph, CodonRestriction
from codeine.graph.nodes import CodonNode, Node
from codeine.graph.tracking import AdvanceResult,BannedSequenceTracker, TrackerState
from codeine.utils.display import format_banned_sequences, format_count, format_restrictions
from codeine.utils.sampling import Sampler, Seedable


NodeState = Tuple[Node, TrackerState]


@dataclass(frozen=True)
class ChoiceResult:
    """
    Cached result of taking one graph choice from one compiled state. The
    "choice" is the graph edge label, i.e. a codon or a context sequence.

    Each ChoiceResult is specific to its location in the graph. The descendant
    counts and weight mass are calculated iteratively by summing the values of
    downstream nodes.
    """
    choice: str
    descendant_count: int
    descendant_weight_mass: float
    next_state: Optional[NodeState]
    is_coding: bool


@dataclass
class CompiledView:
    """
    Cached data for a compiled CodonGraphView, to speed up sampling and enumeration.
    """
    initial_state: NodeState
    n_valid_sequences: int
    choices_by_state: Dict[NodeState, Dict[str, ChoiceResult]]
    codon_pos_by_state: Dict[NodeState, int]
    fixed_choice_by_state: Dict[NodeState, str]
    samplers: dict
    sample_steps: dict


class CodonGraphView:
    """
    View of a codon graph. The view allows optional temporary constraints such as pinned codons
    and banned nucleotide sequences to be added without affecting the underlying codon graph.

    It is on this object that most operations (counting, sampling, enumeration....) take place.
    """

    def __init__(self,
                 graph: CodonGraph,
                 banned_sequences: Optional[Sequence[str]] = None,
                 seed: Seedable = None,
                 rng: Optional[random.Random] = None,
                 ) -> None:
        """
        Constructor for the CodonGraphView

        Parameters
        ----------
        graph
            The underlying codon graph.
        banned_sequences
            Nucleotide sequences that are forbidden in this view.
        seed
            Seed used to initialise a random number generator, if not providing an RNG.
        rng
            Random number generator used by the view for sampling.
        """
        if seed is not None and rng is not None:
            raise ValueError('Provide either seed or rng, not both.')

        self._rng = rng if rng is not None else random.Random(seed)

        self.graph = graph
        self.pinned_codons: Dict[int, List[str]] = {}
        self.banned_sequences: List[str] = self._validate_banned_sequences(banned_sequences)
        self._banned_tracker = None
        self._compiled = None

        self.initial_state = None
        self.n_valid_sequences = None
        self.choices_by_state = {}
        self.codon_pos_by_state = {}
        self.fixed_choice_by_state = {}
        self.samplers = {}
        self.sample_steps = {}

        self.compile()

    @property
    def aa_seq(self):
        """
        The amino acid sequence on the underlying graph.

        Returns
        -------
        The aa seq.
        """
        return self.graph.aa_seq

    def __getitem__(self, index: int) -> str:
        """
        Return the valid sequence at a given index.

        Parameters
        ----------
        index
            Zero-based sequence index.

        Returns
        -------
        str
            The indexed valid DNA sequence.
        """
        return self.sequence_at(index)

    def __iter__(self) -> Generator[str, None, None]:
        """
        Iterate over all valid sequences in this graph view.

        Yields
        ----------
        All valid sequences in the graph, in order.
        """
        yield from self.enumerate()

    def __contains__(self, seq: str) -> bool:
        """
        Does the given seq exist in this space?

        Returns
        ----------
        True if and only if this is a valid sequence in this space.
        """
        return self.contains(seq)

    def __repr__(self) -> str:
        molecule = 'RNA' if self.graph.tt.rna else 'DNA'

        lines = [
            f'{type(self).__name__}',
            '',
            f'Translation table: {self.graph.tt.table_id} ({self.graph.tt.name})',
            f'Molecule type: {molecule}',
            '',
            f'Amino acid sequence ({len(self.aa_seq)} aa)',
            f'{self.aa_seq}',
            ''
        ]

        if self.graph.codon_restrictions:
            lines += [
                'Codon restrictions:',
                *format_restrictions(
                    self.graph.codon_restrictions,
                    label='restricted positions',
                ),
                '',
                ]

        if self.banned_sequences:
            lines += [
                'Banned sequences:',
                *format_banned_sequences(
                    self.banned_sequences,
                ),
                '',
                ]

        if self.pinned_codons:
            lines += [
                'Temporary pins:',
                *format_restrictions(
                    self.pinned_codons,
                    label='pinned positions',
                ),
                '',
                ]

        lines.append(f'Num. valid coding sequences: {format_count(self.n_valid_sequences)}')

        return '\n'.join(lines)

    def pin_codons(self, pinned_codons: Dict[int, CodonRestriction]) -> None:
        """
        Pin (temporarily fix) a codon in this codon graph view

        Parameters
        ----------
        pinned_codons
            A dict specifying which codons to pin, by pos: codon.
        """
        pinned_codons = self.graph.validate_codon_restrictions(pinned_codons)
        self.pinned_codons.update(pinned_codons)
        self.compile()

    def unpin_codons(self, positions: Sequence[int]) -> None:
        """
        Unpin codon nodes by pos.

        Parameters
        ----------
        positions
            A list of positions to unpin.
        """
        for pos in positions:
            if pos < 1 or pos > len(self.graph.codon_nodes):
                raise ValueError(f'Pinned codon position {pos} is out of range.')

            self.pinned_codons.pop(pos, None)

        self.compile()

    def set_pinned_codons(self, pinned_codons: Dict[int, CodonRestriction]) -> None:
        """
        Pin a specified group codons, leaving all others unpinned.

        Parameters
        ----------
        pinned_codons:
            A dict specifying which codons to pin, by pos: codon
        """
        pinned_codons = self.graph.validate_codon_restrictions(pinned_codons)
        self.pinned_codons = dict(pinned_codons)
        self.compile()

    def clear_pins(self) -> None:
        """
        Remove all codon pins from this graph view
        """
        self.pinned_codons.clear()
        self.compile()

    def set_banned_sequences(self, banned_sequences: Sequence[str]) -> None:
        """
        Set banned nucleotide sequences for this view.
        """
        self.banned_sequences = self._validate_banned_sequences(banned_sequences)
        self.compile()

    def contains(self, seq: str) -> bool:
        """
        Check whether a DNA sequence is contained in this view.

        Parameters
        ----------
        seq
            The sequence to check

        Returns
        -------
        True if and only if the sequence is contained in this coding space.
        """
        seq = seq.upper()

        if len(seq) != len(self.graph.aa_seq) * 3:
            return False

        state = self.initial_state
        choices_by_state = self.choices_by_state
        codon_pos_by_state = self.codon_pos_by_state
        fixed_choice_by_state = self.fixed_choice_by_state

        while state is not None:
            pos = codon_pos_by_state.get(state)

            if pos is None:
                choice = fixed_choice_by_state[state]
            else:
                start = (pos - 1) * 3
                choice = seq[start:start + 3]

            result = choices_by_state[state].get(choice)

            if result is None:
                return False

            state = result.next_state

        return True

    def sequence_at(self, index: int) -> str:
        """
        Return the valid sequence at a given index.

        Parameters
        ----------
        index
            Zero-based sequence index.

        Returns
        -------
        str
            The indexed valid DNA sequence.
        """
        if index < 0 or index >= self.n_valid_sequences:
            raise IndexError(f'Sequence index {index} out of range for {self.n_valid_sequences} valid sequences.')

        state = self.initial_state
        sequence = []

        while state is not None:
            result, index = self._choice_result_at_index(state, index)

            if result.is_coding:
                sequence.append(result.choice)

            state = result.next_state

        return ''.join(sequence)

    def sample(self) -> str:
        """
        Sample a DNA sequence from this graph view.
        """
        if self.n_valid_sequences == 0:
            raise ValueError('Cannot sample from an empty coding space.')

        state = self.initial_state
        sequence = []

        while state is not None:
            choice, is_coding, state = self.sample_steps[state].sample()

            if is_coding:
                sequence.append(choice)

        return ''.join(sequence)

    def enumerate(self) -> Generator[str, None, None]:
        """
        Enumerate all valid sequences in this view.

        Yields
        ------
        str
            A valid DNA sequence.
        """
        for index in range(self.n_valid_sequences):
            yield self[index]

    def copy(self) -> 'CodonGraphView':
        """
        Copy this view and all its constraints.

        Returns
        -------
        A copy of the view.
        """
        view = self.graph.view()
        view.pinned_codons = self.pinned_codons.copy()
        view.banned_sequences = self.banned_sequences.copy()
        view.compile()
        return view

    def compile(self) -> None:
        """
        Calculate all graph properties that are derived from its structure plus constraints
        such as pins and banned sequences.

        Remember to do this after editing constraints!
        """
        compiler = ViewCompiler(self)
        compiled = compiler.compile()

        self._banned_tracker = compiler.tracker
        self._compiled = compiled

        self.initial_state = compiled.initial_state
        self.n_valid_sequences = compiled.n_valid_sequences
        self.choices_by_state = compiled.choices_by_state
        self.codon_pos_by_state = compiled.codon_pos_by_state
        self.fixed_choice_by_state = compiled.fixed_choice_by_state
        self.samplers = compiled.samplers
        self.sample_steps = compiled.sample_steps

    def _validate_banned_sequences(self, banned_sequences: Sequence[str]) -> List[str]:
        """
        Check the inputted banned sequences make sense.

        Parameters
        ----------
        banned_sequences
            The list of banned sequences.

        Returns
        -------
        A normalised, de-duplicated list of banned sequences.
        """
        banned_sequences = banned_sequences or []

        normalised = []
        for sequence in banned_sequences:
            sequence = self.graph.normalise_sequence(sequence)

            if len(sequence) == 0:
                raise ValueError('Banned sequences cannot be empty.')

            normalised.append(sequence)

        return sorted(set(normalised))

    def _choices_for_node(self, node) -> List[str]:
        """
        Return choices available to this node in this view.
        """
        if isinstance(node, CodonNode):
            if node.pos in self.pinned_codons:
                return self.pinned_codons[node.pos]

            return node.codons

        return [node.sequence]

    def _choice_result_at_index(self, state: NodeState, index: int) -> Tuple[ChoiceResult, int]:
        """
        Each node + tracking state has a fixed number of descendants and can be enumerated.

        Given a specified index in those descendants, return the choice that contains the
        requested sequence together with the index relative to that chosen subtree.
        """
        results = self.choices_by_state[state]

        if state not in self.codon_pos_by_state:
            result = next(iter(results.values()))
            return result, index

        remaining = index

        for result in results.values():
            if remaining < result.descendant_count:
                return result, remaining

            remaining -= result.descendant_count

        raise RuntimeError('Failed to resolve sequence index.')


class ViewCompiler:
    """
    Compile a CodonGraphView into cached choice, count, and sampling data.
    """

    def __init__(self, view: CodonGraphView) -> None:
        self.view = view
        self.graph = view.graph
        self.tracker = BannedSequenceTracker(view.graph, view.banned_sequences)
        self.advance_cache: Dict[Tuple[Node, TrackerState, str], AdvanceResult] = {}

        self.totals_by_state: Dict[NodeState, Tuple[int, float]] = {}
        self.choices_by_state: Dict[NodeState, Dict[str, ChoiceResult]] = {}
        self.codon_pos_by_state: Dict[NodeState, int] = {}
        self.fixed_choice_by_state: Dict[NodeState, str] = {}

    def compile(self) -> CompiledView:
        """
        Compile descendant counts, graph choices, and samplers.
        """
        initial_state = self._initial_state()
        self._compile_from(initial_state)

        samplers, sample_steps = self._make_samplers()

        return CompiledView(
            initial_state=initial_state,
            n_valid_sequences=self.totals_by_state[initial_state][0],
            choices_by_state=self.choices_by_state,
            codon_pos_by_state=self.codon_pos_by_state,
            fixed_choice_by_state=self.fixed_choice_by_state,
            samplers=samplers,
            sample_steps=sample_steps,
        )

    def _compile_from(self, initial_state: NodeState) -> None:
        """
        Walk the reachable graph states and compile each one after its children.
        """
        initial_node, initial_tracker_state = initial_state
        stack = [(initial_node, initial_tracker_state, False)]

        while stack:
            node, tracker_state, expanded = stack.pop()
            state = self._state(node, tracker_state)

            if state in self.totals_by_state:
                continue

            if node is self.graph.final_node:
                self._compile_final_state(state)
                continue

            if not expanded:
                stack.append((node, tracker_state, True))
                stack.extend(self._uncompiled_children(node, tracker_state))
                continue

            self._compile_state(node, tracker_state)

    def _compile_final_state(self, state: NodeState) -> None:
        """
        Compile the final graph state.
        """
        self.totals_by_state[state] = (1, 1.0)
        self.choices_by_state[state] = {}

    def _compile_state(self, node, tracker_state: TrackerState) -> None:
        """
        Compile one non-final graph state after all valid children have been compiled.
        """
        state = self._state(node, tracker_state)
        choice_results = {}

        descendant_count = 0
        descendant_weight_mass = 0.0

        self._record_state_kind(state, node)

        for choice in self.view._choices_for_node(node):
            result = self._choice_result(node, tracker_state, choice)

            if result is None:
                continue

            choice_results[choice] = result
            descendant_count += result.descendant_count
            descendant_weight_mass += result.descendant_weight_mass

        self.choices_by_state[state] = choice_results
        self.totals_by_state[state] = (descendant_count, descendant_weight_mass)

    def _choice_result(self, node, tracker_state: TrackerState, choice: str) -> Optional[ChoiceResult]:
        """
        Compile the result of taking one outgoing choice from a graph node.
        """
        child = node.transitions.get(choice)

        if child is None:
            return None

        advance = self._advance_tracker(tracker_state, node, choice)

        if advance.banned:
            return None

        child_state = self._state(child, advance.state)
        child_count, child_weight_mass = self.totals_by_state[child_state]

        if child_count == 0:
            return None

        return ChoiceResult(
            choice=choice,
            descendant_count=child_count,
            descendant_weight_mass=self._choice_weight_mass(node, choice, child_weight_mass),
            next_state=None if child is self.graph.final_node else child_state,
            is_coding=isinstance(node, CodonNode),
        )

    def _uncompiled_children(self, node, tracker_state: TrackerState) -> List[Tuple[object, TrackerState, bool]]:
        """
        Return uncompiled child states reachable from a graph state.
        """
        children = []

        for choice in self.view._choices_for_node(node):
            child = node.transitions.get(choice)

            if child is None:
                continue

            advance = self._advance_tracker(tracker_state, node, choice)

            if advance.banned:
                continue

            child_state = self._state(child, advance.state)

            if child_state not in self.totals_by_state:
                children.append((child, advance.state, False))

        return children

    def _make_samplers(self) -> Tuple[dict, dict]:
        """
        Make samplers for each reachable graph state.
        """
        samplers = {}
        sample_steps = {}

        for state, choice_results in self.choices_by_state.items():
            node, _ = state

            if node is self.graph.final_node:
                continue

            runtime_items = []
            runtime_weights = []

            for result in choice_results.values():
                runtime_items.append((result.choice, result.is_coding, result.next_state))
                runtime_weights.append(result.descendant_weight_mass)

            if runtime_items:
                sample_steps[state] = Sampler(runtime_items, runtime_weights, rng=self.view._rng)

            if isinstance(node, CodonNode):
                codons = [result.choice for result in choice_results.values()]
                weights = [result.descendant_weight_mass for result in choice_results.values()]

                if codons:
                    samplers[state] = Sampler(codons, weights, rng=self.view._rng)

        return samplers, sample_steps

    def _record_state_kind(self, state: NodeState, node) -> None:
        """
        Record whether a graph state consumes a codon from the user sequence
        or follows a fixed context sequence.
        """
        if isinstance(node, CodonNode):
            self.codon_pos_by_state[state] = node.pos
        else:
            self.fixed_choice_by_state[state] = node.sequence

    def _choice_weight_mass(self, node, choice: str, child_weight_mass: float) -> float:
        """
        Return the weighted mass for a choice.
        """
        if isinstance(node, CodonNode):
            return self.graph.cw[choice] * child_weight_mass

        return child_weight_mass

    def _initial_state(self) -> NodeState:
        """
        Return the initial compiled graph state.
        """
        return self._state(self.graph.initial_node, self._initial_tracker_state())

    def _initial_tracker_state(self) -> TrackerState:
        """
        Return the initial banned-sequence tracker state.
        """
        if not self.view.banned_sequences:
            return frozenset()

        return self.tracker.initial_state

    def _advance_tracker(self, tracker_state: TrackerState, node: Node, choice: str) -> AdvanceResult:
        """
        Advance banned-sequence tracking after taking a graph step. Results are cached.
        """
        key = (node, tracker_state, choice)

        if key in self.advance_cache:
            return self.advance_cache[key]

        if self.tracker.is_trivial:
            result = AdvanceResult(banned=False, state=tracker_state)
        else:
            step = (node.pos, choice)
            result = self.tracker.advance(step, tracker_state)

        self.advance_cache[key] = result
        return result

    def _state(self, node, tracker_state: TrackerState = frozenset()) -> NodeState:
        """
        Return the compiled state for a graph node plus banned-sequence tracker state.
        """
        return node, tracker_state
