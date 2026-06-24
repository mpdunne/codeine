from dataclasses import dataclass
from typing import Any, Dict, Optional, Tuple


ConstraintState = Any


class PathConstraint:
    """
    Base class for tracking constraints applied while walking a codon graph.
    Designed to track sequence properties that can be calculated by accumulating
    calculations along a path length.

    The idea is to update a state based on the previous state, current node, and choice.
    """

    @property
    def initial_state(self) -> ConstraintState:
        """
        Initial constraint-tracking state.
        """
        return ()

    def advance(
        self,
        state: Any,
        pos: int,
        choice: str,
    ) -> Optional[Any]:
        """
        Advance the constraint state after taking one graph choice.

        Return None if this choice should be rejected.
        """
        return state

    def is_satisfied(self, state: ConstraintState) -> bool:
        """
        Return whether this constraint is satisfied by the current state.
        """
        return True


# nt_diffs, codon_diffs
MutationDistanceState = Tuple[Optional[int], Optional[int]]
MutationDiffCache = Dict[Tuple[int, str], Tuple[int, int]]


@dataclass
class MutationDistanceConstraint(PathConstraint):
    """
    Constrain graph walks by distance from a reference CDS.

    Nucleotide distance counts individual nucleotide changes. Codon distance
    counts codons that differ, regardless of how many nucleotides changed.
    """

    reference_cds: str
    min_nts: Optional[int] = None
    max_nts: Optional[int] = None
    min_codons: Optional[int] = None
    max_codons: Optional[int] = None

    def __post_init__(self) -> None:

        # Store the reference codons once, on init.
        ref_codons = [self.reference_cds[i:i + 3] for i in range(0, len(self.reference_cds), 3)]
        self._ref_codons = tuple(ref_codons)

        # Cache distance calculations for repeated (position, codon) choices.
        self._diff_cache: MutationDiffCache = {}

        self.first_pos = 1
        self.last_pos = len(self._ref_codons)

    @property
    def tracks_nts(self) -> bool:
        """
        Whether nucleotide differences are constrained.
        """
        return self.min_nts is not None or self.max_nts is not None

    @property
    def tracks_codons(self) -> bool:
        """
        Whether codon differences are constrained.
        """
        return self.min_codons is not None or self.max_codons is not None

    @property
    def initial_state(self) -> MutationDistanceState:
        """
        Initial mutation-distance state.
        """
        nt_diffs = 0 if self.tracks_nts else None
        codon_diffs = 0 if self.tracks_codons else None
        return nt_diffs, codon_diffs

    def advance(
        self,
        state: MutationDistanceState,
        pos: int,
        choice: str,
    ) -> Optional[MutationDistanceState]:
        """
        Advance mutation-distance tracking by one graph choice. The state updates on each
        advance by adding the number of nt/codon differences given each next choice.

        Non-codon nodes do not affect distance. Codon nodes add the distance
        between the chosen codon and the reference codon at the same position.
        """
        nt_diffs, codon_diffs = state
        key = (pos, choice)

        if pos < self.first_pos or pos > self.last_pos:
            return state

        cached_diff = self._diff_cache.get(key)
        if cached_diff is None:
            ref_codon = self._ref_codons[pos - 1]
            cached_diff = (
                sum(a != b for a, b in zip(ref_codon, choice)),
                int(ref_codon != choice),
            )
            self._diff_cache[key] = cached_diff

        nt_diff, codon_diff = cached_diff

        if self.tracks_nts:
            nt_diffs += nt_diff

            if self.max_nts is not None and nt_diffs > self.max_nts:
                return None

        if self.tracks_codons:
            codon_diffs += codon_diff

            if self.max_codons is not None and codon_diffs > self.max_codons:
                return None

        return nt_diffs, codon_diffs

    def is_satisfied(self, state: MutationDistanceState) -> bool:
        """
        Check whether a given state satisfies minimum distances.
        """
        nt_diffs, codon_diffs = state

        if self.min_nts is not None and nt_diffs < self.min_nts:
            return False

        if self.min_codons is not None and codon_diffs < self.min_codons:
            return False

        return True
