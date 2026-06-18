from dataclasses import dataclass
from typing import Dict, FrozenSet, List, Sequence, Tuple

from graph.codon import Node


Watch = Tuple[int, int]          # (path_ix, matched_length)
TrackerState = FrozenSet[Watch]


@dataclass(frozen=True)
class BannedPath:
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

    def _find_banned_paths(self) -> Tuple[BannedPath, ...]:
        paths = []

        for sequence in self.banned_sequences:
            for parts, offset in self.graph.find_matching_subpaths(sequence):
                paths.append(
                    BannedPath(
                        sequence=sequence,
                        parts=tuple(parts),
                        offset=offset,
                    )
                )

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
