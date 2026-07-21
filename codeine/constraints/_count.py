import math
from numbers import Real
from typing import Dict, Optional

from codeine.constraints.base import Constraint, ConstraintState, DEAD_STATE, SAFE_STATE
from codeine.graph.base import CodonGraph


INITIAL_STATE = 0


# This constraint is experimental and therefore hidden from the public API
class _CountConstraint(Constraint):
    """
    Base class for constraints that accumulate a non-negative integer count.

    Multiple minimum bounds are combined using the strictest minimum.
    Multiple maximum bounds are combined using the strictest maximum.
    """

    _positions_per_choice: int
    _choice_counts: Dict[str, int]

    def __init__(
        self,
        min_frac: Optional[float] = None,
        max_frac: Optional[float] = None,
        min_perc: Optional[float] = None,
        max_perc: Optional[float] = None,
        min_count: Optional[int] = None,
        max_count: Optional[int] = None,
    ) -> None:
        """
        Parameters
        ----------
        min_frac
            Minimum fraction, between 0 and 1.
        max_frac
            Maximum fraction, between 0 and 1.
        min_perc
            Minimum percentage, between 0 and 100.
        max_perc
            Maximum percentage, between 0 and 100.
        min_count
            Minimum absolute count.
        max_count
            Maximum absolute count.
        """
        self._validate_real("min_frac", min_frac, maximum=1)
        self._validate_real("max_frac", max_frac, maximum=1)
        self._validate_real("min_perc", min_perc, maximum=100)
        self._validate_real("max_perc", max_perc, maximum=100)
        self._validate_count("min_count", min_count)
        self._validate_count("max_count", max_count)

        self.min_frac = min_frac
        self.max_frac = max_frac
        self.min_perc = min_perc
        self.max_perc = max_perc
        self.min_count = min_count
        self.max_count = max_count

        self._n_positions: Optional[int] = None

        self._total_count: Optional[int] = None
        self._resolved_min_count: Optional[int] = None
        self._resolved_max_count: Optional[int] = None

        self._min_remaining = None
        self._max_remaining = None

        self._min_viable_count = None
        self._max_viable_count = None
        self._min_safe_count = None
        self._max_safe_count = None

        self._initial_state = INITIAL_STATE

    @staticmethod
    def _validate_real(name: str, value: Optional[float], maximum: float) -> None:
        """
        Validate a real-valued parameter.
        """
        if value is None:
            return

        if isinstance(value, bool) or not isinstance(value, Real):
            raise TypeError(f'{name} must be a number.')

        if not math.isfinite(value) or not 0 <= value <= maximum:
            raise ValueError(f'{name} must be between 0 and {maximum}.')

    @staticmethod
    def _validate_count(name: str, value: Optional[int]) -> None:
        """
        Validate an integer-valued parameter.
        """
        if value is None:
            return

        if isinstance(value, bool) or not isinstance(value, int):
            raise TypeError(f'{name} must be an integer.')

        if value < 0:
            raise ValueError(f'{name} cannot be negative.')

    @property
    def initial_state(self) -> ConstraintState:
        """
        The initial state of the constraint tracker.
        """
        return self._initial_state

    def link(self, graph: CodonGraph) -> None:
        """
        Link this constraint to a codon graph, resolve count/percentage/fraction
        bounds, and precompute the minimum and maximum remaining counts.
        """
        n_positions = len(graph.codon_nodes)
        self._n_positions = n_positions

        total_count = n_positions * self._positions_per_choice

        self._resolved_min_count = self._resolve_min_count(total_count)
        self._resolved_max_count = self._resolve_max_count(total_count)
        self._total_count = total_count

        # Store the min and max possible remaining counts, in order to shortcut
        # the tracker when we're guaranteed safe / dead.
        self._min_remaining = [0] * (n_positions + 1)
        self._max_remaining = [0] * (n_positions + 1)

        for pos in range(n_positions, 0, -1):
            node = graph.codon_node_by_pos(pos)
            counts = [self._choice_counts[choice] for choice in node.codons]

            ix = pos - 1
            self._min_remaining[ix] = min(counts) + self._min_remaining[ix + 1]
            self._max_remaining[ix] = max(counts) + self._max_remaining[ix + 1]

        # Precompute the accumulated-count thresholds used by advance().
        self._min_viable_count = [0] * (n_positions + 1)
        self._max_viable_count = [0] * (n_positions + 1)
        self._min_safe_count = [0] * (n_positions + 1)
        self._max_safe_count = [0] * (n_positions + 1)

        for pos in range(n_positions + 1):
            self._min_viable_count[pos] = self._resolved_min_count - self._max_remaining[pos]
            self._max_viable_count[pos] = self._resolved_max_count - self._min_remaining[pos]
            self._min_safe_count[pos] = self._resolved_min_count - self._min_remaining[pos]
            self._max_safe_count[pos] = self._resolved_max_count - self._max_remaining[pos]

        min_possible = self._min_remaining[0]
        max_possible = self._max_remaining[0]

        if self._resolved_min_count > self._resolved_max_count:
            self._initial_state = DEAD_STATE
        elif max_possible < self._resolved_min_count or min_possible > self._resolved_max_count:
            self._initial_state = DEAD_STATE
        elif min_possible >= self._resolved_min_count and max_possible <= self._resolved_max_count:
            self._initial_state = SAFE_STATE
        else:
            self._initial_state = INITIAL_STATE

    def advance(
            self,
            state: ConstraintState,
            pos: int,
            choice: str,
    ) -> ConstraintState:
        """
        Return the state after taking a choice at `pos`. Return `DEAD_STATE`
        when the choice violates the constraint. Return `SAFE_STATE` if there are
        no remaining possible sequences that do not pass the constraint.

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

        if state < 0:
            return state

        if pos < 1 or pos > self._n_positions:
            return state

        count = state + self._choice_counts[choice]

        # If there's no way to satisfy the bounds from here, we're dead :(
        if count < self._min_viable_count[pos] or count > self._max_viable_count[pos]:
            return DEAD_STATE

        # If every sequence from here is guaranteed to satisfy the bounds, mark the constraint as safe.
        if self._min_safe_count[pos] <= count <= self._max_safe_count[pos]:
            return SAFE_STATE

        return count

    @property
    def is_trivial(self) -> bool:
        """
        Whether this constraint can never reject any path.
        """
        return self._initial_state == SAFE_STATE

    def _resolve_min_count(self, total_count: int) -> int:
        """
        Resolve the effective minimum count from all specified bounds.
        """
        bounds = [0]

        if self.min_frac is not None:
            bounds.append(math.ceil(self.min_frac * total_count))

        if self.min_perc is not None:
            bounds.append(math.ceil(self.min_perc * total_count / 100))

        if self.min_count is not None:
            bounds.append(self.min_count)

        return max(bounds)

    def _resolve_max_count(self, total_count: int) -> int:
        """
        Resolve the effective maximum count from all specified bounds.
        """
        bounds = [total_count]

        if self.max_frac is not None:
            bounds.append(math.floor(self.max_frac * total_count))

        if self.max_perc is not None:
            bounds.append(math.floor(self.max_perc * total_count / 100))

        if self.max_count is not None:
            bounds.append(self.max_count)

        return min(bounds)
