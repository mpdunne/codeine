from abc import ABC, abstractmethod
from typing import Hashable

from codeine.graph.base import CodonGraph

ConstraintState = Hashable

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
