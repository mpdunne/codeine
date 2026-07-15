from dataclasses import dataclass
from typing import Optional, Tuple

from codeine.constraints.base import Constraint, ConstraintState, DEAD_STATE, SAFE_STATE
from codeine.graph.base import CodonGraph

# nt_diffs, codon_diffs
MutationDistanceState = Tuple[Optional[int], Optional[int]]


@dataclass
class MutationDistanceConstraint(Constraint):
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
        if len(self.reference_cds) % 3:
            raise ValueError('Reference CDS length must be a multiple of three.')

        # Store the reference codons once, on init.
        ref_codons = [self.reference_cds[i:i + 3] for i in range(0, len(self.reference_cds), 3)]
        self._ref_codons = tuple(ref_codons)

        self._tracks_nts = self.min_nts is not None or self.max_nts is not None
        self._tracks_codons = self.min_codons is not None or self.max_codons is not None

        self._initial_state = (
            0 if self._tracks_nts else None,
            0 if self._tracks_codons else None,
        )

        self.first_pos = 1
        self.last_pos = len(self._ref_codons)

    @property
    def tracks_nts(self) -> bool:
        """
        Whether nucleotide differences are constrained.
        """
        return self._tracks_nts

    @property
    def tracks_codons(self) -> bool:
        """
        Whether codon differences are constrained.
        """
        return self._tracks_codons

    @property
    def initial_state(self) -> MutationDistanceState:
        """
        Initial mutation-distance state.
        """
        return self._initial_state

    def advance(
        self,
        state: MutationDistanceState,
        pos: int,
        choice: str,
    ) -> ConstraintState:
        """
        Advance mutation-distance tracking by one graph choice. The state updates on each
        advance by adding the number of nt/codon differences given each next choice.

        Non-codon nodes do not affect distance. Codon nodes add the distance
        between the chosen codon and the reference codon at the same position.
        """
        if state == DEAD_STATE:
            return DEAD_STATE

        if state == SAFE_STATE:
            return SAFE_STATE

        if pos < self.first_pos or pos > self.last_pos:
            return state

        if not self._tracks_nts and not self._tracks_codons:
            return SAFE_STATE

        nt_diffs, codon_diffs = state
        ref_codon = self._ref_codons[pos - 1]

        nt_diff = (
                (ref_codon[0] != choice[0])
                + (ref_codon[1] != choice[1])
                + (ref_codon[2] != choice[2])
        )
        codon_diff = int(nt_diff != 0)

        remaining_codons = self.last_pos - pos

        if self._tracks_nts:
            nt_diffs += nt_diff

            if self.max_nts is not None and nt_diffs > self.max_nts:
                return DEAD_STATE

            if self.min_nts is not None and nt_diffs + 3 * remaining_codons < self.min_nts:
                return DEAD_STATE

            nts_safe = (self.min_nts is None or nt_diffs >= self.min_nts) and \
                       (self.max_nts is None or nt_diffs + 3 * remaining_codons <= self.max_nts)

        else:
            nts_safe = True

        if self._tracks_codons:
            codon_diffs += codon_diff

            if self.max_codons is not None and codon_diffs > self.max_codons:
                return DEAD_STATE

            if self.min_codons is not None and codon_diffs + remaining_codons < self.min_codons:
                return DEAD_STATE

            codons_safe = (self.min_codons is None or codon_diffs >= self.min_codons) and \
                          (self.max_codons is None or codon_diffs + remaining_codons <= self.max_codons)

        else:
            codons_safe = True

        if nts_safe and codons_safe:
            return SAFE_STATE

        return nt_diffs, codon_diffs

    def link(self, graph: CodonGraph) -> None:
        """
        Link up to the graph.
        """
        if len(graph.aa_seq) != len(self._ref_codons):
            raise ValueError('Length of linked graph does not match number of codons for reference CDS.')

    @property
    def is_trivial(self) -> bool:
        """
        Whether this constraint can never reject any path.
        """
        return not self._tracks_nts and not self._tracks_codons
