from typing import Any, Optional

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
