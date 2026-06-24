import math
import random

from dataclasses import dataclass
from typing import Dict, Generator, List, Optional, Sequence, Tuple, Union

from codeine.graph.constraints import ConstraintState, PathConstraint
from codeine.graph.base import CodonGraph, CodonRestriction
from codeine.graph.nodes import CodonNode, Node
from codeine.graph.tracking import AdvanceResult, BannedSequenceTracker, TrackerState
from codeine.utils.display import format_forbidden_motifs, format_count, format_restrictions
from codeine.utils.sampling import Sampler, Seedable


NodeState = Tuple[Node, TrackerState, ConstraintState]


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
    choice_start_by_state: Dict[NodeState, int]
    choice_results_by_state: Dict[NodeState, Tuple[ChoiceResult, ...]]
    codon_pos_by_state: Dict[NodeState, int]
    fixed_choice_by_state: Dict[NodeState, str]
    samplers: dict


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
        self.path_constraint: PathConstraint = PathConstraint()

        self._banned_tracker = BannedSequenceTracker(self.graph, self.banned_sequences)
        self._advance_cache: Dict[Tuple[Node, TrackerState, str], AdvanceResult] = {}

        self._compiled = None
        self._requires_compile = True

        self.initial_state = None
        self.choices_by_state = {}
        self.choice_start_by_state = {}
        self.choice_results_by_state = {}
        self.codon_pos_by_state = {}
        self.fixed_choice_by_state = {}
        self.samplers = {}

    @property
    def aa_seq(self):
        """
        The amino acid sequence on the underlying graph.

        Returns
        -------
        The aa seq.
        """
        return self.graph.aa_seq

    @property
    def translation_table(self):
        """
        The translation table used by the codon graph.
        """
        return self.graph.tt

    @property
    def codon_weights(self):
        """
        The codon weights used by the codon graph.
        """
        return self.graph.cw

    @property
    def codon_restrictions(self):
        """
        Fixed codon restrictions.
        """
        return self.graph.codon_restrictions

    @property
    def context_l(self):
        """
        The left context sequence.
        """
        return self.graph.context_l

    @property
    def context_r(self):
        """
        The right context sequence.
        """
        return self.graph.context_r

    @property
    def n_valid_sequences(self) -> int:
        """
        Number of valid coding sequences in this view.
        """
        if self._requires_compile:
            self.compile()

        return self._compiled.n_valid_sequences

    def __getitem__(self, index: Union[int, slice]) -> Union[str, List[str]]:
        """
        Return one valid sequence, or a list of valid sequences for a slice.

        Parameters
        ----------
        index
            Zero-based sequence index, or slice of sequence indices.

        Returns
        -------
        str or list of str
            The indexed valid DNA sequence, or a list of valid DNA sequences.
        """
        if isinstance(index, slice):
            return self.sequences_at(index)

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
        if self._requires_compile:
            self.compile()

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
                *format_forbidden_motifs(
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

        lines.append(f'Num. valid coding sequences: {format_count(self._compiled.n_valid_sequences)}')

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
        self._requires_compile = True

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

        self._requires_compile = True

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
        self._requires_compile = True

    def clear_pins(self) -> None:
        """
        Remove all codon pins from this graph view
        """
        self.pinned_codons.clear()
        self._requires_compile = True

    def set_banned_sequences(self, banned_sequences: Sequence[str]) -> None:
        """
        Set banned nucleotide sequences for this view.

        Banned-sequence tracking depends only on the graph and banned sequences,
        not on temporary pins, so it is rebuilt only when the banned list changes.
        """
        self.banned_sequences = self._validate_banned_sequences(banned_sequences)
        self._banned_tracker = BannedSequenceTracker(self.graph, self.banned_sequences)
        self._advance_cache.clear()
        self._requires_compile = True

    def set_path_constraint(self, path_constraint: PathConstraint) -> None:
        """
        Set an additional generic path constraint for this view.

        The view does not interpret the constraint. It only lets the constraint
        track state while walking the graph and reject choices or final states.
        """
        self.path_constraint = path_constraint
        self._requires_compile = True

    def clear_path_constraint(self) -> None:
        """
        Remove the additional generic path constraint from this view.
        """
        self.path_constraint = PathConstraint()
        self._requires_compile = True

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
        if self._requires_compile:
            self.compile()

        seq = seq.upper()

        if len(seq) != len(self.graph.aa_seq) * 3:
            return False

        state = self.initial_state
        choices_by_state = self.choices_by_state
        choice_start_by_state = self.choice_start_by_state
        fixed_choice_by_state = self.fixed_choice_by_state

        while state is not None:
            start = choice_start_by_state.get(state)

            if start is None:
                choice = fixed_choice_by_state[state]
            else:
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
        if self._requires_compile:
            self.compile()

        if index < 0 or index >= self._compiled.n_valid_sequences:
            raise IndexError(
                f'Sequence index {index} out of range for {self._compiled.n_valid_sequences} valid sequences.'
            )

        state = self.initial_state
        sequence = []

        choice_results_by_state = self.choice_results_by_state
        codon_pos_by_state = self.codon_pos_by_state

        while state is not None:
            results = choice_results_by_state[state]

            if state not in codon_pos_by_state:
                result = results[0]
            else:
                remaining = index

                for result in results:
                    if remaining < result.descendant_count:
                        index = remaining
                        break

                    remaining -= result.descendant_count
                else:
                    raise RuntimeError('Failed to resolve sequence index.')

                sequence.append(result.choice)

            state = result.next_state

        return ''.join(sequence)

    def sample(self) -> str:
        """
        Sample a DNA sequence from this graph view.
        """
        if self._requires_compile:
            self.compile()

        if self._compiled.n_valid_sequences == 0:
            raise ValueError('Cannot sample from an empty coding space.')

        state = self.initial_state
        sequence = []
        samplers = self.samplers

        while state is not None:
            choice, is_coding, state = samplers[state].sample()

            if is_coding:
                sequence.append(choice)

        return ''.join(sequence)

    def sequences_at(self, index_slice: slice) -> List[str]:
        """
        Return valid sequences from a slice.

        Parameters
        ----------
        index_slice
            Slice of zero-based sequence indices.

        Returns
        -------
        list of str
            The sliced valid DNA sequences.
        """
        if self._requires_compile:
            self.compile()

        start, stop, step = index_slice.indices(self._compiled.n_valid_sequences)

        if step != 1:
            return [self.sequence_at(index) for index in range(start, stop, step)]

        if start == 0:
            sequences = []
            for index, sequence in enumerate(self.enumerate()):
                if index >= stop:
                    break

                sequences.append(sequence)

            return sequences

        return [*self.enumerate_range(start, stop)]

    def enumerate(self) -> Generator[str, None, None]:
        """
        Enumerate all valid sequences in this view.

        Yields
        ------
        str
            A valid DNA sequence.
        """
        if self._requires_compile:
            self.compile()

        stack = [(self.initial_state, [])]

        while stack:
            state, sequence_parts = stack.pop()

            if state is None:
                yield ''.join(sequence_parts)
                continue

            results = self.choice_results_by_state[state]

            if not results:
                continue

            if state not in self.codon_pos_by_state:
                stack.append((results[0].next_state, sequence_parts))
                continue

            for result in reversed(results):
                stack.append((
                    result.next_state,
                    [*sequence_parts, result.choice],
                ))

    def enumerate_range(self, start: int = 0, stop: Optional[int] = None) -> Generator[str, None, None]:
        """
        Enumerate valid sequences from start up to, but not including, stop.
        """
        if self._requires_compile:
            self.compile()

        n_sequences = self._compiled.n_valid_sequences

        if stop is None:
            stop = n_sequences

        if start < 0 or stop < start or stop > n_sequences:
            raise IndexError('Enumeration range is out of bounds.')

        for index in range(start, stop):
            yield self.sequence_at(index)

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
        view.path_constraint = self.path_constraint
        view._banned_tracker = self._banned_tracker
        view._advance_cache = self._advance_cache.copy()

        view._compiled = self._compiled
        view._requires_compile = self._requires_compile

        view.initial_state = self.initial_state
        view.choices_by_state = self.choices_by_state
        view.choice_start_by_state = self.choice_start_by_state
        view.choice_results_by_state = self.choice_results_by_state
        view.codon_pos_by_state = self.codon_pos_by_state
        view.fixed_choice_by_state = self.fixed_choice_by_state
        view.samplers = self.samplers

        return view

    def compile(self) -> None:
        """
        Calculate all graph properties that are derived from its structure plus constraints
        such as pins and banned sequences.

        Remember to do this after editing constraints!
        """
        compiler = ViewCompiler(self)
        compiled = compiler.compile()

        self._compiled = compiled
        self._requires_compile = False

        self.initial_state = compiled.initial_state
        self.choices_by_state = compiled.choices_by_state
        self.choice_results_by_state = compiled.choice_results_by_state
        self.choice_start_by_state = compiled.choice_start_by_state
        self.codon_pos_by_state = compiled.codon_pos_by_state
        self.fixed_choice_by_state = compiled.fixed_choice_by_state
        self.samplers = compiled.samplers

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
            sequence = self.translation_table.normalise_sequence(sequence)

            if len(sequence) == 0:
                raise ValueError('Banned sequences cannot be empty.')

            normalised.append(sequence)

        return sorted(set(normalised))


class ViewCompiler:
    """
    Compile a CodonGraphView into cached choice, count, and sampling data.
    """

    def __init__(self, view: CodonGraphView) -> None:
        self.view = view
        self.graph = view.graph
        self.tracker = view._banned_tracker
        self.advance_cache = view._advance_cache
        self.path_constraint = view.path_constraint

        self.totals_by_state: Dict[NodeState, Tuple[int, float]] = {}
        self.choices_by_state: Dict[NodeState, Dict[str, ChoiceResult]] = {}
        self.choice_start_by_state: Dict[NodeState, int] = {}
        self.choice_results_by_state: Dict[NodeState, Tuple[ChoiceResult, ...]] = {}
        self.codon_pos_by_state: Dict[NodeState, int] = {}
        self.fixed_choice_by_state: Dict[NodeState, str] = {}

        self.choices_by_node = {
            node: tuple(self._get_choices_for_node(node))
            for node in self.graph.nodes
            if node is not self.graph.final_node
        }

    def compile(self) -> CompiledView:
        """
        Compile descendant counts, graph choices, and samplers.
        """
        initial_state = self._initial_state()
        self._compile_from(initial_state)

        self._compile_choice_result_tuples()
        self._compile_choice_starts()
        samplers = self._make_samplers()

        return CompiledView(
            initial_state=initial_state,
            n_valid_sequences=self.totals_by_state[initial_state][0],
            choices_by_state=self.choices_by_state,
            choice_results_by_state=self.choice_results_by_state,
            choice_start_by_state=self.choice_start_by_state,
            codon_pos_by_state=self.codon_pos_by_state,
            fixed_choice_by_state=self.fixed_choice_by_state,
            samplers=samplers,
        )

    def _compile_from(self, initial_state: NodeState) -> None:
        """
        Walk the reachable graph states and compile each one after its children.
        """
        initial_node, initial_tracker_state, initial_constraint_state = initial_state
        stack = [(initial_node, initial_tracker_state, initial_constraint_state, False)]

        while stack:
            node, tracker_state, constraint_state, expanded = stack.pop()
            state = self._state(node, tracker_state, constraint_state)

            if state in self.totals_by_state:
                continue

            if node is self.graph.final_node:
                self._compile_final_state(state, constraint_state)
                continue

            if not expanded:
                stack.append((node, tracker_state, constraint_state, True))
                stack.extend(self._uncompiled_children(node, tracker_state, constraint_state))
                continue

            self._compile_state(node, tracker_state, constraint_state)

    def _compile_final_state(self, state: NodeState, constraint_state: ConstraintState) -> None:
        """
        Compile the final graph state.
        """
        if self.path_constraint.is_satisfied(constraint_state):
            self.totals_by_state[state] = (1, 0.0)
        else:
            self.totals_by_state[state] = (0, -math.inf)

        self.choices_by_state[state] = {}

    def _compile_state(
        self,
        node,
        tracker_state: TrackerState,
        constraint_state: ConstraintState,
    ) -> None:
        """
        Compile one non-final graph state after all valid children have been compiled.
        """
        state = self._state(node, tracker_state, constraint_state)
        raw_results = []

        self._record_state_kind(state, node)

        for choice in self.choices_by_node[node]:
            result = self._choice_result(node, tracker_state, constraint_state, choice)

            if result is None:
                continue

            raw_results.append(result)

        results = raw_results

        choice_results = {}
        descendant_count = 0
        descendant_weight_masses = []

        for result in results:
            choice_results[result.choice] = result
            descendant_count += result.descendant_count
            descendant_weight_masses.append(result.descendant_weight_mass)

        descendant_weight_mass = self._sum_weight_masses(descendant_weight_masses)

        self.choices_by_state[state] = choice_results
        self.totals_by_state[state] = (descendant_count, descendant_weight_mass)

    def _compile_choice_result_tuples(self) -> None:
        """
        Store each state's choice results as tuples for fast indexed traversal.
        """
        self.choice_results_by_state = {
            state: tuple(choice_results.values())
            for state, choice_results in self.choices_by_state.items()
        }

    def _get_choices_for_node(self, node) -> List[str]:
        """
        Return choices available to this node in this view.
        """
        if isinstance(node, CodonNode):
            if node.pos in self.view.pinned_codons:
                return self.view.pinned_codons[node.pos]

            return node.codons

        return [node.sequence]

    def _choice_result(
        self,
        node,
        tracker_state: TrackerState,
        constraint_state: ConstraintState,
        choice: str,
    ) -> Optional[ChoiceResult]:
        """
        Compile the result of taking one outgoing choice from a graph node.
        """
        child = node.transitions.get(choice)

        if child is None:
            return None

        advance = self._advance_tracker(tracker_state, node, choice)

        if advance.banned:
            return None

        next_constraint_state = self.path_constraint.advance(
            constraint_state,
            node.pos,
            choice,
        )

        if next_constraint_state is None:
            return None

        child_state = self._state(child, advance.state, next_constraint_state)
        child_count, child_weight_mass = self.totals_by_state[child_state]

        if child_count == 0:
            return None

        descendant_weight_mass = self._choice_weight_mass(node, choice, child_weight_mass)

        if descendant_weight_mass == -math.inf:
            return None

        return ChoiceResult(
            choice=choice,
            descendant_count=child_count,
            descendant_weight_mass=descendant_weight_mass,
            next_state=None if child is self.graph.final_node else child_state,
            is_coding=isinstance(node, CodonNode),
        )

    def _uncompiled_children(
        self,
        node,
        tracker_state: TrackerState,
        constraint_state: ConstraintState,
    ) -> List[Tuple[object, TrackerState, ConstraintState, bool]]:
        """
        Return uncompiled child states reachable from a graph state.
        """
        children = []

        for choice in self.choices_by_node[node]:
            child = node.transitions.get(choice)

            if child is None:
                continue

            advance = self._advance_tracker(tracker_state, node, choice)

            if advance.banned:
                continue

            next_constraint_state = self.path_constraint.advance(
                constraint_state,
                node.pos,
                choice,
            )

            if next_constraint_state is None:
                continue

            child_state = self._state(child, advance.state, next_constraint_state)

            if child_state not in self.totals_by_state:
                children.append((child, advance.state, next_constraint_state, False))

        return children

    def _normalise_choice_results(self, results: List[ChoiceResult]) -> List[ChoiceResult]:
        """
        Rescale descendant weight masses within one state.

        Only relative weights matter for sampling, so this prevents long paths
        from underflowing toward zero.
        """
        max_mass = max((result.descendant_weight_mass for result in results), default=0.0)

        if max_mass <= 0:
            return results

        return [
            ChoiceResult(
                choice=result.choice,
                descendant_count=result.descendant_count,
                descendant_weight_mass=result.descendant_weight_mass / max_mass,
                next_state=result.next_state,
                is_coding=result.is_coding,
            )
            for result in results
        ]

    def _normalise_weights(self, weights: List[float]) -> List[float]:
        """
        Rescale sampler weights defensively.
        """
        if not weights:
            return weights

        max_weight = max(weights)

        if max_weight == -math.inf:
            return [1.0] * len(weights)

        return [math.exp(weight - max_weight) for weight in weights]

    def _make_samplers(self) -> dict:
        """
        Make samplers for each reachable graph state.
        """
        samplers = {}

        for state, choice_results in self.choice_results_by_state.items():
            node, _, _ = state

            if node is self.graph.final_node:
                continue

            runtime_items = []
            runtime_weights = []

            for result in choice_results:
                runtime_items.append((result.choice, result.is_coding, result.next_state))
                runtime_weights.append(result.descendant_weight_mass)

            if runtime_items:
                runtime_weights = self._normalise_weights(runtime_weights)
                samplers[state] = Sampler(runtime_items, runtime_weights, rng=self.view._rng)

        return samplers

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
            weight = self.graph.cw[choice]

            if weight <= 0:
                return -math.inf

            return math.log(weight) + child_weight_mass

        return child_weight_mass

    def _sum_weight_masses(self, weights: List[float]) -> float:
        """
        Sum weight masses defensively.
        """
        weights = [weight for weight in weights if weight != -math.inf]

        if not weights:
            return -math.inf

        max_weight = max(weights)

        return max_weight + math.log(
            sum(math.exp(weight - max_weight) for weight in weights)
        )

    def _initial_state(self) -> NodeState:
        """
        Return the initial compiled graph state.
        """
        return self._state(
            self.graph.initial_node,
            self._initial_tracker_state(),
            self.path_constraint.initial_state,
        )

    def _initial_tracker_state(self) -> TrackerState:
        """
        Return the initial banned-sequence tracker state.
        """
        if not self.view.banned_sequences:
            return frozenset()

        return self.tracker.initial_state

    def _compile_choice_starts(self) -> None:
        """
        Store sequence slice starts for coding states.
        """
        self.choice_start_by_state = {
            state: (pos - 1) * 3
            for state, pos in self.codon_pos_by_state.items()
        }

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

    def _state(
            self,
            node,
            tracker_state: TrackerState = frozenset(),
            constraint_state: ConstraintState = (),
    ) -> NodeState:
        """
        Return the compiled state for a graph node plus tracker states.
        """
        return node, tracker_state, constraint_state
