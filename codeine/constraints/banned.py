from typing import Dict, FrozenSet, List, NamedTuple, Optional, Sequence, Tuple

from codeine.graph.base import CodonGraph
from codeine.graph.nodes import CodonNode

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

# The tracker state is a set of watches. We update the
# watches every time we make a choice.
BannedTrackerState = FrozenSet[Watch]

# Internal transition value:
#   None   -> banned sequence completed
#   Watch  -> continue watching this path
TransitionValue = Optional[Watch]

# Integerised version of a banned state
BannedTrackerStateId = int


class AdvanceResult(NamedTuple):
    """
    The result of moving the tracker state forward, indicating whether we are
    currently in a disallowed state, and if not, what the new state is.
    """
    banned: bool
    state: BannedTrackerState = frozenset()


class AdvanceIdResult(NamedTuple):
    """
    The result of moving the tracker state forward, indicating whether we are
    currently in a disallowed state, and if not, what the new state is.
    """
    banned: bool
    state_id: BannedTrackerStateId = 0


CLEAR_ADVANCE_RESULT = AdvanceResult(banned=False)
BANNED_ADVANCE_RESULT = AdvanceResult(banned=True)

CLEAR_ADVANCE_ID_RESULT = AdvanceIdResult(banned=False, state_id=0)
BANNED_ADVANCE_ID_RESULT = AdvanceIdResult(banned=True, state_id=0)


class BannedSequenceTracker:
    """
    Tracks progress along concrete "banned" graph subpaths.

    A SubPath stores:

        sequence
            The sequence being tracked.

        steps
            Concrete graph emissions as (pos, choice) pairs.

        offset
            Where `sequence` starts inside the first choice.

    State is a frozenset of watches:

        (path_ix, matched_length)

    meaning that `matched_length` bases of paths[path_ix].sequence
    have already matched.

    Transitions are precomputed:

        choice -> (path_ix, matched_length) -> banned | next watch | dead
    """

    def __init__(self, graph: CodonGraph, banned_sequences: Sequence[str]) -> None:
        """
        Constructor for the BannedSequenceTracker class.

        Parameters
        ----------
        graph
            The codon graph on which to operate.
        banned_sequences
            A collection of "banned" sequences that we need to watch out for.
        """
        self.graph = graph
        self.banned_sequences = tuple(sequence.upper() for sequence in banned_sequences)
        self.initial_state: BannedTrackerState = frozenset()

        self.state_ids: Dict[BannedTrackerState, BannedTrackerStateId] = {self.initial_state: 0}
        self.states: List[BannedTrackerState] = [self.initial_state]

        self.advance_id_cache: Dict[Tuple[Step, BannedTrackerStateId], AdvanceIdResult] = {}

        self.paths = self._find_banned_paths()
        self.starts = self._build_starts()
        self.transitions = self._build_transitions()

    @property
    def is_trivial(self) -> bool:
        """
        Whether this tracker is trivial - i.e. there are no paths that would generate
        a sequence containing a banned sequence.

        Returns
        -------
        True if and only if the tracker is trivial.
        """
        return len(self.paths) == 0

    def _find_banned_paths(self) -> Tuple[SubPath, ...]:
        """
        Find every concrete graph path that can generate a banned sequence.

        Each returned SubPath records the emitted sequence, the graph steps
        required to produce it, and the offset at which the banned sequence begins
        within the first emitted choice.

        Returns
        -------
        Tuple[SubPath, ...]
            All graph subpaths capable of producing one of the banned sequences.
        """
        paths = []

        for sequence in self.banned_sequences:
            paths.extend(_find_matching_subpaths(self.graph, sequence))

        return tuple(paths)

    def _build_starts(self) -> Dict[Step, Tuple[TransitionValue, ...]]:
        """
        Build the initial watch transitions for each possible graph step.

        The returned mapping records which watches should be created when a given
        step is taken. If a banned sequence is completed immediately, the transition
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
        Build transitions between tracker states.

        For each emitted graph choice, records how every active watch should
        advance. A transition value of None indicates that the banned sequence
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

    def advance(
            self,
            step: Step,
            state: BannedTrackerState,
    ) -> AdvanceResult:
        """
        Move the tracker state forward after taking a graph step.

        Parameters
        ----------
        step
            The graph step just taken, as (graph pos, choice).
        state
            The current tracker state.

        Returns
        -------
        AdvanceResult
            Whether the step completed a banned sequence, and otherwise the
            updated tracker state.
        """
        starts = self.starts.get(step)

        if starts is None and not state:
            return CLEAR_ADVANCE_RESULT

        _pos, choice = step
        transitions = self.transitions.get(choice)
        next_state = set()

        if starts is not None:
            for watch in starts:
                if watch is None:
                    return BANNED_ADVANCE_RESULT

                next_state.add(watch)

        if transitions is not None:
            for watch in state:
                if watch not in transitions:
                    continue

                next_watch = transitions[watch]

                if next_watch is None:
                    return BANNED_ADVANCE_RESULT

                next_state.add(next_watch)

        if not next_state:
            return CLEAR_ADVANCE_RESULT

        return AdvanceResult(banned=False, state=frozenset(next_state))

    def advance_id(
            self,
            step: Step,
            state_id: BannedTrackerStateId,
    ) -> AdvanceIdResult:
        key = (step, state_id)

        cached = self.advance_id_cache.get(key)
        if cached is not None:
            return cached

        state = self.states[state_id]
        starts = self.starts.get(step)

        if starts is None and not state:
            self.advance_id_cache[key] = CLEAR_ADVANCE_ID_RESULT
            return CLEAR_ADVANCE_ID_RESULT

        _pos, choice = step
        transitions = self.transitions.get(choice)
        next_state = set()

        if starts is not None:
            for watch in starts:
                if watch is None:
                    self.advance_id_cache[key] = BANNED_ADVANCE_ID_RESULT
                    return BANNED_ADVANCE_ID_RESULT

                next_state.add(watch)

        if transitions is not None:
            for watch in state:
                next_watch = transitions.get(watch)

                if next_watch is None:
                    if watch in transitions:
                        self.advance_id_cache[key] = BANNED_ADVANCE_ID_RESULT
                        return BANNED_ADVANCE_ID_RESULT

                    continue

                next_state.add(next_watch)

        if not next_state:
            self.advance_id_cache[key] = CLEAR_ADVANCE_ID_RESULT
            return CLEAR_ADVANCE_ID_RESULT

        result = AdvanceIdResult(
            banned=False,
            state_id=self._get_or_register_state_id(frozenset(next_state)),
        )

        self.advance_id_cache[key] = result
        return result

    def _get_or_register_state_id(
            self,
            state: BannedTrackerState,
    ) -> BannedTrackerStateId:
        state_id = self.state_ids.get(state)

        if state_id is None:
            state_id = len(self.states)
            self.state_ids[state] = state_id
            self.states.append(state)

        return state_id


def _find_matching_subpaths(
        graph: CodonGraph,
        sequence: str,
) -> List[SubPath]:
    """
    Find all graph subpaths that can emit a given banned sequence.

    Parameters
    ----------
    graph
        The codon graph to search.
    sequence
        The banned sequence to search for.

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

        for choice, child in node.transitions.items():
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
