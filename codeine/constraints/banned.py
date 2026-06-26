from dataclasses import dataclass
from typing import Dict, FrozenSet, List, Sequence, Tuple

from codeine.graph.base import CodonGraph
from codeine.graph.nodes import Node, CodonNode

# A step is a decision in the codon graph, i.e. (graph pos, choice)
Step = Tuple[int, str]


@dataclass(frozen=True)
class SubPath:
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


@dataclass(frozen=True)
class AdvanceResult:
    """
    The result of moving the tracker state forward, indicating whether we are
    currently in a disallowed state, and if not, what the new state is.
    """
    banned: bool
    state: BannedTrackerState = frozenset()


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

        (path_ix, matched_length), choice -> banned | next watch | dead
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

        self.paths = self._find_banned_paths()
        self.starts = self._build_starts()
        self.transitions = self._build_transitions()

    @property
    def is_trivial(self) -> bool:
        return len(self.paths) == 0

    def _find_banned_paths(self) -> Tuple[SubPath, ...]:
        paths = []

        for sequence in self.banned_sequences:
            for parts, offset in _find_matching_subpaths(self.graph, sequence):
                steps = tuple(
                    (node.pos, choice)
                    for node, choice in parts
                )

                paths.append(
                    SubPath(sequence=sequence, steps=steps, offset=offset)
                )

        return tuple(paths)

    def _build_starts(self) -> Dict[Step, Tuple[AdvanceResult, ...]]:
        starts: Dict[Step, List[AdvanceResult]] = {}

        for path_ix, path in enumerate(self.paths):
            first_step = path.steps[0]
            _pos, first_choice = first_step

            emitted = first_choice[path.offset:]
            matched_length = min(len(emitted), len(path.sequence))

            if matched_length >= len(path.sequence):
                result = AdvanceResult(banned=True)
            else:
                state = frozenset({(path_ix, matched_length)})
                result = AdvanceResult(banned=False, state=state)

            starts.setdefault(first_step, []).append(result)

        return {
            key: tuple(results)
            for key, results in starts.items()
        }

    def _build_transitions(self) -> Dict[Tuple[Watch, str], AdvanceResult]:
        transitions = {}

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
                    transitions[(watch, choice)] = AdvanceResult(banned=True)
                    break

                if remaining.startswith(choice):
                    matched_length += len(choice)
                    state = frozenset({(path_ix, matched_length)})
                    transitions[(watch, choice)] = AdvanceResult(banned=False, state=state)
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
        _pos, choice = step
        next_state = set()

        for result in self.starts.get(step, ()):
            if result.banned:
                return AdvanceResult(banned=True)

            next_state.update(result.state)

        for watch in state:
            result = self.transitions.get((watch, choice))

            if result is None:
                continue

            if result.banned:
                return AdvanceResult(banned=True)

            next_state.update(result.state)

        return AdvanceResult(banned=False, state=frozenset(next_state))


MatchedPath = Tuple[Tuple[Node, str], ...]
MatchedSubPath = Tuple[MatchedPath, int]


def _find_matching_subpaths(graph: CodonGraph, sequence: str) \
        -> List[MatchedSubPath]:
    """
    For a given sequence, find subpaths in the graph that match that sequence.
    Return each found subpath in the following format:

        (
            [
                (node1, choice_1),
                (node2, choice_2),
                ...
            ]
            offset,  # Where the path starts relative to first node's codon choice.
        )

    Parameters
    ----------
    sequence
        The sequence to search for.

    Returns
    -------
    A tuple consisting of a list of (node, sequence) pairs, and the start offset for matching the sequence.
    """

    sequence = sequence.upper()

    if len(sequence) == 0:
        raise ValueError('Sequence cannot be empty.')

    matches: List[MatchedSubPath] = []
    candidate_matches = []

    # First, check which nodes we can start at.
    for node in graph.nodes:
        if node is graph.end_node:
            continue

        for choice, child in node.transitions.items():
            for offset in range(len(choice)):
                choice_subsequence = choice[offset:]

                if choice_subsequence.startswith(sequence):
                    # Bingo!
                    matches.append((((node, choice),), offset))

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
                matches.append((partial_path, offset))
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
                        matches.append((
                            partial_path + ((node, choice),),
                            offset,
                        ))

                    else:
                        # Hard luck this time.
                        continue

        return reinspect

    while candidate_matches:
        candidate_matches = reinspect_candidate_matches(candidate_matches)

    return matches