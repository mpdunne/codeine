from abc import ABC, abstractmethod
from typing import Dict, List, Tuple, NamedTuple, FrozenSet, Optional

from codeine.constraints.base import Constraint, ConstraintState, DEAD_STATE, SAFE_STATE
from codeine.graph.base import CodonGraph


# A step is a decision in the codon graph, i.e. (graph pos, choice)
Step = Tuple[int, str]


class SubPath(NamedTuple):
    """
    A subpath in the codon graph, indicating a sequence that can be obtained
    by following a specified sequence of steps, starting at a given offset.
    """
    sequence: str
    steps: Tuple[Step, ...]
    offset: int


# A watch (path_ix, matched_length) is the status of a single
# watched path, indicating how much of the path has been seen so far.
Watch = Tuple[int, int]

# The constraint state is a set of watches. We update the
# watches every time we make a choice.
TrackerState = FrozenSet[Watch]

# Integer ID for a registered tracker state.
TrackerStateId = int

# Internal transition value:
#   None   -> forbidden subpath completed
#   Watch  -> continue watching this path
TransitionValue = Optional[Watch]


class SubPathConstraint(Constraint, ABC):
    """
    Base class for constraints that forbid specific graph subpaths.

    Subclasses identify the forbidden subpaths for a linked graph. This class
    tracks partial matches against those paths and rejects any sequence that
    completes one of them.
    """

    def __init__(self) -> None:

        initial_tracker_state: TrackerState = frozenset()
        self.initial_state_id: int = 0

        self.state_ids: Dict[TrackerState, TrackerStateId] = {initial_tracker_state: self.initial_state_id}
        self.states: List[TrackerState] = [initial_tracker_state]

        self.advance_cache: Dict[Tuple[Step, TrackerStateId], ConstraintState] = {}

        self.graph: Optional[CodonGraph] = None
        self.aa_seq = None
        self.context_l = None
        self.context_r = None
        self.translation_table = None

        self.paths: Tuple[SubPath, ...] = ()
        self.starts: Dict[Step, Tuple[TransitionValue, ...]] = {}
        self.transitions: Dict[str, Dict[Watch, TransitionValue]] = {}

    @property
    def initial_state(self) -> ConstraintState:
        """
        Return the registered empty tracker-state ID.
        """
        return self.initial_state_id

    def advance(
            self,
            state: ConstraintState,
            pos: int,
            choice: str,
    ) -> ConstraintState:
        """
        Advance the constraint state after taking one graph choice.

        Returns DEAD_STATE if the choice completes a forbidden subpath; otherwise
        returns the integer ID of the updated tracker state.
        """
        if state == DEAD_STATE:
            return DEAD_STATE

        if state == SAFE_STATE:
            return SAFE_STATE

        state_id = state
        step = (pos, choice)
        key = (step, state_id)

        cached = self.advance_cache.get(key)
        if cached is not None:
            return cached

        tracker_state = self.states[state_id]
        starts = self.starts.get(step)

        if starts is None and not tracker_state:
            self.advance_cache[key] = self.initial_state_id
            return self.initial_state_id

        transitions = self.transitions.get(choice)
        next_state = set()

        if starts is not None:
            for watch in starts:
                if watch is None:
                    self.advance_cache[key] = DEAD_STATE
                    return DEAD_STATE

                next_state.add(watch)

        if transitions is not None:
            for watch in tracker_state:
                next_watch = transitions.get(watch)

                if next_watch is None:
                    if watch in transitions:
                        self.advance_cache[key] = DEAD_STATE
                        return DEAD_STATE

                    continue

                next_state.add(next_watch)

        if not next_state:
            self.advance_cache[key] = self.initial_state_id
            return self.initial_state_id

        next_state_id = self._get_or_register_state_id(frozenset(next_state))
        self.advance_cache[key] = next_state_id
        return next_state_id

    def link(self, graph: CodonGraph) -> None:
        """
        Link the constraint to a graph.
        """
        self._reset_states()

        self.graph = graph
        self.aa_seq = graph.aa_seq
        self.context_l = graph.context_l
        self.context_r = graph.context_r
        self.translation_table = graph.tt

        self.paths = self._find_paths()
        self.starts = self._build_starts()
        self.transitions = self._build_transitions()

    @property
    def is_trivial(self) -> bool:
        """
        Whether this constraint is trivial - i.e. there are no paths that would
        generate a sequence containing a forbidden subpath.

        Returns
        -------
        True if and only if the constraint is trivial.
        """
        return len(self.paths) == 0

    def _reset_states(self) -> None:
        """
        Reset registered tracker states and cached transitions.
        """
        initial_state: TrackerState = frozenset()

        self.state_ids = {initial_state: self.initial_state_id}
        self.states = [initial_state]
        self.advance_cache.clear()

    @abstractmethod
    def _find_paths(self) -> Tuple[SubPath, ...]:
        """
        Find the forbidden subpaths for the linked graph.

        Returns
        -------
        Tuple[SubPath, ...]
            The concrete graph subpaths that must not be completed.
        """
        raise NotImplementedError

    def _build_starts(self) -> Dict[Step, Tuple[TransitionValue, ...]]:
        """
        Build the initial watch transitions for each possible graph step.

        The returned mapping records which watches should be created when a given
        step is taken. If a forbidden motif is completed immediately, the transition
        value is None.

        Returns
        -------
        Dict[Step, Tuple[TransitionValue, ...]]
            Mapping from graph step to the watches that should be started after
            taking that step.
        """
        starts: Dict[Step, List[TransitionValue]] = {}

        for path_ix, path in enumerate(self.paths):
            first_step = path.steps[0]
            _pos, first_choice = first_step

            emitted = first_choice[path.offset:]
            matched_length = min(len(emitted), len(path.sequence))

            if matched_length >= len(path.sequence):
                result = None
            else:
                result = (path_ix, matched_length)

            starts.setdefault(first_step, []).append(result)

        return {
            key: tuple(results)
            for key, results in starts.items()
        }

    def _build_transitions(self) -> Dict[str, Dict[Watch, TransitionValue]]:
        """
        Build transitions between active watches.

        For each emitted graph choice, records how every active watch should
        advance. A transition value of None indicates that the forbidden motif
        has been completed.

        Returns
        -------
        Dict[str, Dict[Watch, TransitionValue]]
            Mapping from emitted graph choice to the corresponding watch
            transitions.
        """
        transitions: Dict[str, Dict[Watch, TransitionValue]] = {}

        for path_ix, path in enumerate(self.paths):
            matched_length = min(
                len(path.steps[0][1]) - path.offset,
                len(path.sequence),
            )

            if matched_length >= len(path.sequence):
                continue

            for _pos, choice in path.steps[1:]:
                watch = (path_ix, matched_length)
                remaining = path.sequence[matched_length:]

                if choice.startswith(remaining):
                    transitions.setdefault(choice, {})[watch] = None
                    break

                if remaining.startswith(choice):
                    matched_length += len(choice)
                    transitions.setdefault(choice, {})[watch] = (path_ix, matched_length)
                    continue

                break

        return transitions

    def _get_or_register_state_id(
            self,
            state: TrackerState,
    ) -> TrackerStateId:
        """
        Return the integer ID for a tracker state, creating one if needed.

        Parameters
        ----------
        state
            The concrete frozenset-of-watches tracker state.

        Returns
        -------
        TrackerStateId
            Stable integer ID for the given tracker state.
        """
        state_id = self.state_ids.get(state)

        if state_id is None:
            state_id = len(self.states)
            self.state_ids[state] = state_id
            self.states.append(state)

        return state_id
