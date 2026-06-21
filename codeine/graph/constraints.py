from typing import Any, Optional

from codeine.graph.nodes import Node


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