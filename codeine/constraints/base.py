from abc import ABC, abstractmethod
from typing import Hashable, Optional, Tuple

from codeine.graph.base import CodonGraph
from codeine.graph.factors import Factor

ConstraintState = Hashable

SAFE_STATE = -2
DEAD_STATE = -1


class Constraint(ABC):
    """
    Base class for tracking constraints applied while walking a codon graph.
    Designed to track sequence properties that can be calculated by accumulating
    calculations along a path length.

    The idea is to update a state based on the previous state, current node, and choice.
    """

    @property
    @abstractmethod
    def initial_state(self) -> ConstraintState:
        """
        Initial constraint-tracking state.
        """
        pass

    @abstractmethod
    def advance(
        self,
        state: ConstraintState,
        pos: int,
        choice: str,
    ) -> ConstraintState:
        """
        Return the state after taking a choice at `pos`.
        Return `DEAD_STATE` when the choice violates the constraint.

        Parameters
        ----------
        state
            The input state.
        pos
            The current position.
        choice
            The codon choice.

        Returns
        -------
        An updated state.
        """
        pass

    @abstractmethod
    def link(self, graph: CodonGraph) -> None:
        """
        Link up this constraint with a codon graph and precompute any relevant data.

        Parameters
        ----------
        graph
            The graph to link.
        """
        pass

    def factors(self) -> Optional[Tuple[Factor, ...]]:
        """
        Describe the constraint as a collection of local rules.

        Each factor checks only the graph positions involved in one rule. A
        sequence satisfies this constraint if and only if every factor accepts
        the sequence's choices. This exposes the dependencies needed by
        a factorised compiler to solve independent groups of rules separately.

        Call link before requesting factors. The existing flat compiler continues
        to use initial_state and advance; implementing this method is optional.

        Returns
        -------
        Tuple of factors, or None if this representation is unsupported.
        An empty tuple means the constraint permits every graph assignment.
        """
        return None

    @property
    def is_trivial(self) -> bool:
        """
        Whether this constraint can never reject any path.
        """
        return False
