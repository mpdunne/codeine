import math
import random

from bisect import bisect_right
from typing import Mapping, Optional

from codeine.graph.factor_count import ComponentModelCounter


class FactorSampler:
    """
    Sample valid assignments using a reusable AND/OR plan.

    The counter supplies the decomposition and cached transitions. This class
    adds probability masses and a sampling plan, keeping weights separate from
    validity counts. Independent groups are sampled separately; within a group,
    each choice is weighted by the mass of all its valid completions.

    The plan contains no random generator. Callers can share it while retaining
    their own random streams, and samplers with different weights can share a
    counter without affecting one another.
    """

    def __init__(
        self,
        counter: ComponentModelCounter,
        weights: Optional[Mapping[int, Mapping[str, float]]] = None,
    ) -> None:
        """
        Parameters
        ----------
        counter
            The factor model to sample from.
        weights
            Choice weights by graph position. Missing weights default to 1.
            Weights must be finite and non-negative; zero excludes a choice
            from sampling without removing it from the valid sequence count.
        """
        self.counter = counter
        weights = weights or {}
        self._log_weights = {}

        for pos, choices in counter.domains.items():
            log_weights = {}

            for choice in choices:
                weight = float(weights.get(pos, {}).get(choice, 1.0))

                if not math.isfinite(weight) or weight < 0:
                    raise ValueError('Choice weights must be finite and non-negative.')

                log_weights[choice] = math.log(weight) if weight > 0 else -math.inf

            self._log_weights[pos] = log_weights

        self._free_log_masses = {
            pos: self._logsumexp(log_weights.values())
            for pos, log_weights in self._log_weights.items()
        }
        self._free_samplers = {
            pos: self._make_sampler(tuple(log_weights), tuple(log_weights.values()))
            for pos, log_weights in self._log_weights.items()
        }
        self._mass_cache = {}
        self._plan_cache = {}
        self._sampling_plan = None

    @staticmethod
    def _logsumexp(values):
        """
        Add probability masses in log space without overflowing their products.
        """
        values = tuple(values)
        maximum = max(values, default=-math.inf)

        if maximum == -math.inf:
            return -math.inf

        return maximum + math.log(sum(math.exp(value - maximum) for value in values))

    @staticmethod
    def _make_sampler(items, log_masses):
        """
        Cache positive choices and cumulative weights scaled by the largest mass.
        """
        maximum = max(log_masses, default=-math.inf)
        kept_items = []
        cumulative = []
        total = 0.0

        for item, log_mass in zip(items, log_masses):
            if log_mass == -math.inf:
                continue

            weight = math.exp(log_mass - maximum)

            if weight == 0:
                continue

            total += weight
            kept_items.append(item)
            cumulative.append(total)

        return tuple(kept_items), tuple(cumulative), total

    @staticmethod
    def _draw(sampler, rng):
        items, cumulative, total = sampler

        if len(items) == 1:
            return items[0]

        index = bisect_right(cumulative, rng.random() * total)
        return items[min(index, len(items) - 1)]

    def _free_log_mass(self, positions):
        return sum(self._free_log_masses[pos] for pos in self.counter._positions(positions))

    def _mass_component(self, component_id):
        counter = self.counter
        stack = [(component_id, False)]

        while stack:
            current, expanded = stack.pop()

            if current in self._mass_cache:
                continue

            roots = counter._component_roots[current]

            if counter.manager.FALSE in roots:
                result = -math.inf
            elif not roots:
                result = 0.0
            else:
                components = counter._split_components(current)

                if len(components) > 1:
                    children = tuple((child, 0.0) for child in components)
                else:
                    pos = counter._choose_variable(current)
                    children = tuple(
                        (child, self._log_weights[pos][choice] + self._free_log_mass(free))
                        for choice in counter.domains[pos]
                        for child, free in (counter._advance_component(current, pos, choice),)
                    )

                if not expanded:
                    stack.append((current, True))
                    stack.extend((child, False) for child, _weight in children)
                    continue

                if len(components) > 1:
                    result = sum(self._mass_cache[child] for child, _weight in children)
                else:
                    result = self._logsumexp(self._mass_cache[child] + weight for child, weight in children)

            self._mass_cache[current] = result

        return self._mass_cache[component_id]

    def _plan_component(self, component_id):
        counter = self.counter
        stack = [(component_id, False)]

        while stack:
            current, expanded = stack.pop()

            if current in self._plan_cache:
                continue

            roots = counter._component_roots[current]

            if not roots:
                self._plan_cache[current] = ('and', ())
                continue

            components = counter._split_components(current)

            if len(components) > 1:
                if not expanded:
                    stack.append((current, True))
                    stack.extend((child, False) for child in components)
                    continue

                plan = ('and', tuple(components))
            else:
                pos = counter._choose_variable(current)
                branches = []
                log_masses = []

                for choice in counter.domains[pos]:
                    child, free = counter._advance_component(current, pos, choice)
                    log_mass = self._log_weights[pos][choice] + self._free_log_mass(free) + self._mass_component(child)

                    if log_mass == -math.inf:
                        continue

                    branches.append((choice, child, tuple(counter._positions(free))))
                    log_masses.append(log_mass)

                if not expanded:
                    stack.append((current, True))
                    stack.extend((child, False) for _choice, child, _free in branches)
                    continue

                items = tuple(branches)
                plan = ('or', pos, self._make_sampler(items, log_masses))

            self._plan_cache[current] = plan

        return component_id

    def prepare_sampling(self):
        """
        Build the plan once, rejecting empty spaces and zero total weight.
        """
        if self._sampling_plan is None:
            counter = self.counter

            if counter.count() == 0:
                raise ValueError('Cannot sample from an empty factor model.')

            log_mass = (
                self._mass_component(counter._initial_component_id)
                + self._free_log_mass(counter._initial_free_variables)
            )

            if log_mass == -math.inf:
                raise ValueError('Cannot sample from a factor model with zero total weight.')

            self._sampling_plan = (
                self._plan_component(counter._initial_component_id),
                tuple(counter._positions(counter._initial_free_variables)),
            )

        return self._sampling_plan

    def sample(self, rng: random.Random, n: Optional[int] = None):
        """
        Sample one assignment, or a batch, using the caller's random generator.

        Assignments map graph positions to chosen codons or context sequences.
        """
        plan = self.prepare_sampling()

        if n is None:
            return self._sample(plan, rng)

        if n < 0:
            raise ValueError('n must be non-negative.')

        return [self._sample(plan, rng) for _ in range(n)]

    def _sample(self, plan, rng):
        assignment = {}
        stack = [plan]

        while stack:
            node_id, free_positions = stack.pop()
            node = self._plan_cache[node_id]

            for pos in free_positions:
                assignment[pos] = self._draw(self._free_samplers[pos], rng)

            if node[0] == 'and':
                stack.extend((child, ()) for child in reversed(node[1]))
            else:
                _kind, pos, sampler = node
                choice, child, free_positions = self._draw(sampler, rng)
                assignment[pos] = choice
                stack.append((child, free_positions))

        return assignment
