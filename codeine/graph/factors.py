from typing import Dict, FrozenSet, Mapping, Protocol, Sequence, Tuple


# A relation is (reference pos, compare pos, matching choices).
# For example, (1, 2, {'AAA': frozenset({'AAG', 'AAA'})}) matches when
# position 1 is AAA and position 2 is either AAG or AAA.
ChoiceRelation = Tuple[int, int, Dict[str, FrozenSet[str]]]


class Factor(Protocol):
    """
    A rule exposing its dependent positions and how to build its decision diagram.
    """
    scope: FrozenSet[int]

    def compile(self, diagrams) -> int:
        """
        Return the root of a diagram that accepts assignments satisfying this rule.
        """
        pass


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

    def compile(self, diagrams) -> int:
        """
        Build a diagram accepting every choice combination this factor permits.
        """
        return diagrams.allowed_factor(self.relations)

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


class CountFactor:
    """
    Require the sum of per-position choice counts to lie within inclusive bounds.

    This describes GC counts or mutation distances without exposing their meaning
    to the compiler. Unlike a forbidden-choice factor, it accepts a bounded sum.
    """

    def __init__(self, counts: Mapping[int, Mapping[str, int]], minimum: int, maximum: int) -> None:
        self.counts = {pos: dict(choices) for pos, choices in counts.items()}
        self.minimum = minimum
        self.maximum = maximum
        self.scope = frozenset(counts)

    def compile(self, diagrams) -> int:
        """
        Build the bounded-sum decision diagram using the available graph choices.
        """
        return diagrams.bounded_sum(self.counts, self.minimum, self.maximum)
