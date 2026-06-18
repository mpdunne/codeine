from dataclasses import dataclass
from typing import Dict, FrozenSet, List, Sequence, Tuple

from codeine.graph.graph import Node, CodonNode, CodonGraph


Watch = Tuple[int, int]          # (path_ix, matched_length)
TrackerState = FrozenSet[Watch]


@dataclass(frozen=True)
class SubPath:
    sequence: str
    parts: Tuple[Tuple[Node, str], ...]
    offset: int


@dataclass(frozen=True)
class AdvanceResult:
    banned: bool
    state: TrackerState = frozenset()


class BannedSequenceTracker:
    """
    Tracks progress along concrete banned graph paths.

    State is a frozenset of watches:

        (path_ix, matched_length)

    meaning that `matched_length` bases of paths[path_ix].sequence
    have already matched on the current graph walk.
    """

    def __init__(self, graph, banned_sequences: Sequence[str]) -> None:
        self.graph = graph
        self.banned_sequences = tuple(banned_sequences)
        self.initial_state: TrackerState = frozenset()

        self.paths = self._find_banned_paths()
        self.starts = self._build_starts()

    @property
    def is_trivial(self) -> bool:
        return len(self.paths) == 0

    def _find_banned_paths(self) -> Tuple[SubPath, ...]:
        paths = []

        for sequence in self.banned_sequences:
            for parts, offset in _find_matching_subpaths(self.graph, sequence):
                subpath = SubPath(sequence=sequence, parts=tuple(parts), offset=offset)
                paths.append(subpath)

        return tuple(paths)

    def _build_starts(self) -> Dict[Tuple[Node, str], Tuple[Watch, ...]]:
        starts: Dict[Tuple[Node, str], List[Watch]] = {}

        for path_ix, path in enumerate(self.paths):
            first_node, first_choice = path.parts[0]

            emitted = first_choice[path.offset:]
            matched_length = min(len(emitted), len(path.sequence))

            starts.setdefault((first_node, first_choice), []).append(
                (path_ix, matched_length)
            )

        return {
            key: tuple(watches)
            for key, watches in starts.items()
        }

    def advance(
        self,
        node: Node,
        state: TrackerState,
        choice: str,
    ) -> AdvanceResult:
        """
        Advance after taking `choice` from `node`.

        Returns an AdvanceResult describing whether a banned
        sequence has been completed and which watches remain active.
        """
        next_state = set()

        # Start newly possible watches.
        for path_ix, matched_length in self.starts.get((node, choice), ()):
            path = self.paths[path_ix]

            if matched_length >= len(path.sequence):
                return AdvanceResult(banned=True)

            next_state.add((path_ix, matched_length))

        # Continue existing watches.
        for path_ix, matched_length in state:
            path = self.paths[path_ix]
            remaining = path.sequence[matched_length:]

            if not remaining.startswith(choice):
                continue

            new_matched_length = matched_length + len(choice)

            if new_matched_length >= len(path.sequence):
                return AdvanceResult(banned=True)

            next_state.add((path_ix, new_matched_length))

        return AdvanceResult(
            banned=False,
            state=frozenset(next_state),
        )


def _find_matching_subpaths(graph: CodonGraph, sequence: str) \
        -> List[Tuple[List[Tuple[Node, str]], int]]:
    """
    For a given sequence, find subpaths in the graph that match that sequence.
    Return each found subpath in the following format:

        (
            [
                (node1, codon_1),
                (node2, codon_2),
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

    matches = []
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
                    matches.append(([(node, choice)], offset))

                elif sequence.startswith(choice_subsequence):
                    # Maybe bingo! Maygo!
                    candidate_matches.append(([(node, choice)], offset, len(choice_subsequence)))

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

            if node is graph.final_node:
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
                        reinspect.append((partial_path + [(node, choice)], offset, seen_length + choice_length))

                    else:
                        # Hard luck this time.
                        continue

                else:

                    if choice.startswith(remaining_sequence):
                        # Wahoo!
                        matches.append((partial_path + [(node, choice)], offset))

                    else:
                        # Hard luck this time.
                        continue

        return reinspect

    while candidate_matches:
        candidate_matches = reinspect_candidate_matches(candidate_matches)

    return matches
