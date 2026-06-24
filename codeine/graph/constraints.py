from typing import Any, Optional, Tuple
from dataclasses import dataclass

from codeine.graph.nodes import CodonNode, Node


ConstraintState = Any


class PathConstraint:
    """
    Optional constraint applied while walking a codon graph.

    This is deliberately generic: the graph view does not need to know what the
    constraint means. It only asks whether a graph choice is allowed and whether
    the final state is acceptable.
    """

    @property
    def initial_state(self) -> ConstraintState:
        """
        Initial constraint-tracking state.
        """
        return ()

    def advance(
        self,
        state: ConstraintState,
        node: Node,
        choice: str,
    ) -> Optional[ConstraintState]:
        """
        Advance the constraint state after taking one graph choice.

        Return None if this choice should be rejected.
        """
        return state

    def accepts_final(self, state: ConstraintState) -> bool:
        """
        Return whether a completed graph walk is accepted.
        """
        return True


# nt_diffs, codon_diffs
MutationDistanceState = Tuple[int, int]


@dataclass()
class MutationDistanceConstraint(PathConstraint):
    """
    Constrain graph walks by nucleotide and/or codon distance from a reference CDS.
    """

    reference_cds: str
    min_nts: Optional[int] = None
    max_nts: Optional[int] = None
    min_codons: Optional[int] = None
    max_codons: Optional[int] = None

    def __post_init__(self) -> None:
        ref_codons = [self.reference_cds[i:i + 3] for i in range(0, len(self.reference_cds), 3)]
        self._ref_codons = tuple(ref_codons)
        self._diff_cache = {}

    @property
    def tracks_nts(self) -> bool:
        """
        Whether nucleotide differences should be tracked.
        """
        return self.min_nts is not None or self.max_nts is not None

    @property
    def tracks_codons(self) -> bool:
        """
        Whether codon differences should be tracked.
        """
        return self.min_codons is not None or self.max_codons is not None

    @property
    def initial_state(self) -> MutationDistanceState:
        """
        Initial mutation-distance state.

        Counts start at zero for each distance type being tracked.
        Distance types that are not constrained are stored as None.
        """
        nt_diffs = 0 if self.tracks_nts else None
        codon_diffs = 0 if self.tracks_codons else None
        return nt_diffs, codon_diffs

    def advance(
        self,
        state: ConstraintState,
        node,
        choice: str,
    ) -> Optional[MutationDistanceState]:
        """
        Advance mutation-distance tracking by one graph step.

        When a codon node is traversed, nucleotide and codon differences
        relative to the reference CDS are accumulated.

        Returns
        -------
        MutationDistanceState
            Updated mutation-distance state.

        None
            If a maximum-distance constraint has been exceeded.
        """
        if not isinstance(node, CodonNode):
            return state

        nt_diffs, codon_diffs = state
        key = (node.pos, choice)

        cached = self._diff_cache.get(key)
        if cached is None:
            ref_codon = self._ref_codons[node.pos - 1]
            cached = (
                sum(a != b for a, b in zip(ref_codon, choice)),
                int(ref_codon != choice),
            )
            self._diff_cache[key] = cached

        nt_diff, codon_diff = cached

        if self.tracks_nts:
            nt_diffs += nt_diff

            if self.max_nts is not None and nt_diffs > self.max_nts:
                return None

        if self.tracks_codons:
            codon_diffs += codon_diff

            if self.max_codons is not None and codon_diffs > self.max_codons:
                return None

        return nt_diffs, codon_diffs

    def accepts_final(self, state: ConstraintState) -> bool:
        """
        Check whether a completed sequence satisfies the minimum
        mutation-distance constraints.
        """
        nt_diffs, codon_diffs = state

        if self.min_nts is not None and nt_diffs < self.min_nts:
            return False

        if self.min_codons is not None and codon_diffs < self.min_codons:
            return False

        return True
