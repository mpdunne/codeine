from typing import Dict, FrozenSet, List, NamedTuple, Optional, Sequence, Tuple, Union

from codeine.constraints.base import Constraint, ConstraintState, DEAD_STATE, SAFE_STATE
from codeine.graph.base import CodonGraph
from codeine.graph.nodes import CodonNode
from codeine.motifs.restriction import RestrictionSite

# A step is a decision in the codon graph, i.e. (graph pos, choice)
Step = Tuple[int, str]


class SubPath(NamedTuple):
    """
    A subpath in the codon graph, indicating a sequence that can be obtained
    by following a specified sequence of steps, starting at a given offset.
    """
    sequence: str
    steps: Tuple[Step, ...]
    offset: int


# A watch (path_ix, matched_length) is the status of a single
# watched path, indicating how much of the path has been seen so far.
Watch = Tuple[int, int]

# The constraint state is a set of watches. We update the
# watches every time we make a choice.
TrackerState = FrozenSet[Watch]

# Integer ID for a registered tracker state.
TrackerStateId = int

# Internal transition value:
#   None   -> forbidden motif completed
#   Watch  -> continue watching this path
TransitionValue = Optional[Watch]


Motif = Union[str, RestrictionSite]
Motifs = Union[Motif, Sequence[Motif]]


class ForbiddenMotifConstraint(Constraint):
    """
    A constraint class for preventing specified nucleotide motifs from
    occurring in graph sequences.

    Motifs may be provided as nucleotide strings or as ``RestrictionSite``
    objects. Restriction sites are expanded into their concrete recognition
    sequences during construction.
    """

    def __init__(self, forbidden_motifs: Motifs) -> None:
        """
        Parameters
        ----------
        forbidden_motifs
            One or more nucleotide motifs or restriction sites that must not occur.
        """
        if isinstance(forbidden_motifs, (str, RestrictionSite)):
            forbidden_motifs = [forbidden_motifs]

        sequences = []

        for motif in forbidden_motifs:
            if isinstance(motif, RestrictionSite):
                sequences.extend(motif.motifs)
                continue

            if not isinstance(motif, str):
                raise TypeError('Forbidden motifs must be strings or RestrictionSite objects.')

            motif = motif.upper()

            if motif == '':
                raise ValueError('Forbidden motifs cannot be empty.')

            if not set(motif) <= set('ACGTU'):
                raise ValueError('Forbidden motifs must contain only A, C, G, T, or U.')

            sequences.append(motif)

        self.forbidden_sequences = tuple(sorted(set(sequences)))

        initial_tracker_state: TrackerState = frozenset()
        self.initial_state_id: int = 0

        self.state_ids: Dict[TrackerState, TrackerStateId] = {initial_tracker_state: self.initial_state_id}
        self.states: List[TrackerState] = [initial_tracker_state]

        self.advance_cache: Dict[Tuple[Step, TrackerStateId], ConstraintState] = {}

        self.graph: Optional[CodonGraph] = None
        self.paths: Tuple[SubPath, ...] = ()
        self.starts: Dict[Step, Tuple[TransitionValue, ...]] = {}
        self.transitions: Dict[str, Dict[Watch, TransitionValue]] = {}

    @property
    def initial_state(self) -> ConstraintState:
        """
        Return the registered empty tracker-state ID.
        """
        return self.initial_state_id

    @property
    def is_trivial(self) -> bool:
        """
        Whether this constraint is trivial - i.e. there are no paths that would
        generate a sequence containing a forbidden motif.

        Returns
        -------
        True if and only if the constraint is trivial.
        """
        return len(self.paths) == 0

    def link(self, graph: CodonGraph) -> None:
        """
        Link the constraint to a graph.
        """
        initial_state: TrackerState = frozenset()

        self.state_ids = {initial_state: self.initial_state_id}
        self.states = [initial_state]
        self.advance_cache.clear()

        self.forbidden_sequences = tuple(sorted({
            graph.tt.normalise_sequence(sequence)
            for sequence in self.forbidden_sequences
        }))

        self.graph = graph
        self.paths = self._find_matching_paths()
        self.starts = self._build_starts()
        self.transitions = self._build_transitions()

    def advance(
            self,
            state: ConstraintState,
            pos: int,
            choice: str,
    ) -> ConstraintState:
        """
        Advance the constraint state after taking one graph choice.

        Returns DEAD_STATE if the choice completes a forbidden motif; otherwise
        returns the integer ID of the updated tracker state.
        """
        if state == DEAD_STATE:
            return DEAD_STATE

        if state == SAFE_STATE:
            return SAFE_STATE

        state_id = state
        step = (pos, choice)
        key = (step, state_id)

        cached = self.advance_cache.get(key)
        if cached is not None:
            return cached

        tracker_state = self.states[state_id]
        starts = self.starts.get(step)

        if starts is None and not tracker_state:
            self.advance_cache[key] = self.initial_state_id
            return self.initial_state_id

        transitions = self.transitions.get(choice)
        next_state = set()

        if starts is not None:
            for watch in starts:
                if watch is None:
                    self.advance_cache[key] = DEAD_STATE
                    return DEAD_STATE

                next_state.add(watch)

        if transitions is not None:
            for watch in tracker_state:
                next_watch = transitions.get(watch)

                if next_watch is None:
                    if watch in transitions:
                        self.advance_cache[key] = DEAD_STATE
                        return DEAD_STATE

                    continue

                next_state.add(next_watch)

        if not next_state:
            self.advance_cache[key] = self.initial_state_id
            return self.initial_state_id

        next_state_id = self._get_or_register_state_id(frozenset(next_state))
        self.advance_cache[key] = next_state_id
        return next_state_id

    def _find_matching_paths(self) -> Tuple[SubPath, ...]:
        """
        Find every concrete graph path that can generate a forbidden motif.

        Each returned SubPath records the emitted sequence, the graph steps
        required to produce it, and the offset at which the forbidden motif begins
        within the first emitted choice.

        Returns
        -------
        Tuple[SubPath, ...]
            All graph subpaths capable of producing one of the forbidden motifs.
        """
        paths = []

        for sequence in self.forbidden_sequences:
            paths.extend(_find_matching_subpaths(self.graph, sequence))

        return tuple(paths)

    def _build_starts(self) -> Dict[Step, Tuple[TransitionValue, ...]]:
        """
        Build the initial watch transitions for each possible graph step.

        The returned mapping records which watches should be created when a given
        step is taken. If a forbidden motif is completed immediately, the transition
        value is None.

        Returns
        -------
        Dict[Step, Tuple[TransitionValue, ...]]
            Mapping from graph step to the watches that should be started after
            taking that step.
        """
        starts: Dict[Step, List[TransitionValue]] = {}

        for path_ix, path in enumerate(self.paths):
            first_step = path.steps[0]
            _pos, first_choice = first_step

            emitted = first_choice[path.offset:]
            matched_length = min(len(emitted), len(path.sequence))

            if matched_length >= len(path.sequence):
                result = None
            else:
                result = (path_ix, matched_length)

            starts.setdefault(first_step, []).append(result)

        return {
            key: tuple(results)
            for key, results in starts.items()
        }

    def _build_transitions(self) -> Dict[str, Dict[Watch, TransitionValue]]:
        """
        Build transitions between active watches.

        For each emitted graph choice, records how every active watch should
        advance. A transition value of None indicates that the forbidden motif
        has been completed.

        Returns
        -------
        Dict[str, Dict[Watch, TransitionValue]]
            Mapping from emitted graph choice to the corresponding watch
            transitions.
        """
        transitions: Dict[str, Dict[Watch, TransitionValue]] = {}

        for path_ix, path in enumerate(self.paths):
            matched_length = min(
                len(path.steps[0][1]) - path.offset,
                len(path.sequence),
            )

            if matched_length >= len(path.sequence):
                continue

            for _pos, choice in path.steps[1:]:
                watch = (path_ix, matched_length)
                remaining = path.sequence[matched_length:]

                if choice.startswith(remaining):
                    transitions.setdefault(choice, {})[watch] = None
                    break

                if remaining.startswith(choice):
                    matched_length += len(choice)
                    transitions.setdefault(choice, {})[watch] = (path_ix, matched_length)
                    continue

                break

        return transitions

    def _get_or_register_state_id(
            self,
            state: TrackerState,
    ) -> TrackerStateId:
        """
        Return the integer ID for a tracker state, creating one if needed.

        Parameters
        ----------
        state
            The concrete frozenset-of-watches tracker state.

        Returns
        -------
        TrackerStateId
            Stable integer ID for the given tracker state.
        """
        state_id = self.state_ids.get(state)

        if state_id is None:
            state_id = len(self.states)
            self.state_ids[state] = state_id
            self.states.append(state)

        return state_id


###############################################################################
# Algorithm for finding graph subpaths matching a nucleotide sequence.
###############################################################################


def _find_matching_subpaths(
        graph: CodonGraph,
        sequence: str,
) -> List[SubPath]:
    """
    Find all graph subpaths that can emit a given sequence.

    Parameters
    ----------
    graph
        The codon graph to search.
    sequence
        The sequence to search for.

    Returns
    -------
    Matching subpaths as a list of SubPath objects.
    """
    sequence = sequence.upper()

    if len(sequence) == 0:
        raise ValueError('Sequence cannot be empty.')

    matches: List[SubPath] = []
    candidate_matches = []

    def add_match(partial_path, offset):
        steps = tuple(
            (node.pos, choice)
            for node, choice in partial_path
        )

        matches.append(
            SubPath(
                sequence=sequence,
                steps=steps,
                offset=offset,
            )
        )

    # First, check which nodes we can start at.
    for node in graph.nodes:
        if node is graph.end_node:
            continue

        for choice in node.transitions:
            for offset in range(len(choice)):
                choice_subsequence = choice[offset:]

                if choice_subsequence.startswith(sequence):
                    # Bingo!
                    add_match(((node, choice),), offset)

                elif sequence.startswith(choice_subsequence):
                    # Maybe bingo! Maygo!
                    candidate_matches.append((
                        ((node, choice),),
                        offset,
                        len(choice_subsequence),
                    ))

    def reinspect_candidate_matches(candidate_matches):
        reinspect = []

        for partial_path, offset, seen_length in candidate_matches:
            previous_node, previous_choice = partial_path[-1]
            node = previous_node.transitions[previous_choice]

            remaining_sequence = sequence[seen_length:]

            if remaining_sequence == '':
                # Fantastic!
                add_match(partial_path, offset)
                continue

            if node is graph.end_node:
                continue

            if isinstance(node, CodonNode):
                choice_length = 3
            else:
                choice_length = len(next(iter(node.transitions)))

            # Sneaky shortcut if we've crossed into the right context:
            if isinstance(node, CodonNode):
                pos = node.pos

                remaining_sequence_length = len(sequence) - seen_length
                remaining_coding_length = 3 * (len(graph.aa_seq) - pos + 1)

                if remaining_sequence_length > remaining_coding_length:
                    sequence_end = sequence[seen_length + remaining_coding_length:]

                    if not graph.context_r.startswith(sequence_end):
                        continue

            for choice, child in node.transitions.items():

                if len(remaining_sequence) >= choice_length:

                    if remaining_sequence.startswith(choice):
                        # Keep going...
                        reinspect.append((
                            partial_path + ((node, choice),),
                            offset,
                            seen_length + choice_length,
                        ))

                    else:
                        # Hard luck this time.
                        continue

                else:

                    if choice.startswith(remaining_sequence):
                        # Wahoo!
                        add_match(
                            partial_path + ((node, choice),),
                            offset,
                        )

                    else:
                        # Hard luck this time.
                        continue

        return reinspect

    while candidate_matches:
        candidate_matches = reinspect_candidate_matches(candidate_matches)

    return matches
