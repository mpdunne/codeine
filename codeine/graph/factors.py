from typing import Dict, FrozenSet, Mapping, Sequence, Tuple


# A relation is (reference pos, compare pos, matching choices).
# For example, (1, 2, {'AAA': frozenset({'AAG', 'AAA'})}) matches when
# position 1 is AAA and position 2 is either AAG or AAA.
ChoiceRelation = Tuple[int, int, Dict[str, FrozenSet[str]]]


class ChoiceFactor:
    """
    One forbidden combination of graph choices, expressed as matching relations.

    A constraint can produce several factors. For example, each possible repeat
    location produces one factor: the repeat occurs only if all the nucleotide
    comparisons at that location match. Any one factor rejecting a sequence is
    enough to violate the constraint.

    The scope is the set of graph positions used by this factor. Exposing these
    positions and their matching choices lets a compiler find independent groups
    of rules without knowing whether they describe repeats, motifs, or another
    constraint. Factors may overlap; they are not necessarily independent.
    """

    def __init__(self, relations: Sequence[ChoiceRelation]) -> None:
        """
        Parameters
        ----------
        relations
            Tuples of (reference position, compare position, matching choices).
            Each mapping gives the compare choices that match a reference choice.
            The factor rejects only when every relation matches. Positions are
            graph node positions, including context nodes.
        """
        self.relations = tuple(relations)

        self.scope = frozenset(
            pos
            for reference_pos, compare_pos, _choices in self.relations
            for pos in (reference_pos, compare_pos)
        )

    def rejects(self, assignments: Mapping[int, str]) -> bool:
        """
        Return whether the chosen graph steps complete this forbidden combination.

        An empty set of relations rejects every assignment. A choice missing
        from a relation's mapping does not match that relation.

        Parameters
        ----------
        assignments
            Mapping from graph position to chosen codon or context sequence.
            Must include every position in this factor's scope.
        """
        return all(
            assignments[compare_pos] in choices.get(assignments[reference_pos], ())
            for reference_pos, compare_pos, choices in self.relations
        )
