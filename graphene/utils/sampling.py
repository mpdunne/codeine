import bisect
import random

from typing import Any, List, Hashable, Union


class Sampler:
    """
    A precomputed sampler to speed up weighted sampling.
    """
    def __init__(self, items: List[Any], weights: List[Union[int, float]] = None, seed: Hashable = None):
        if not items:
            raise ValueError("Items cannot be empty")

        if weights is None:
            weights = [1] * len(items)

        if len(items) != len(weights):
            raise ValueError("Items and weights must have same length")

        total = sum(weights)
        if total <= 0:
            raise ValueError("Weights must sum to > 0")

        self.items = list(items)
        self.rng = random.Random(seed)

        cumulative = []
        running = 0

        for weight in weights:
            if weight < 0:
                raise ValueError("weights cannot be negative")

            running += weight
            cumulative.append(running / total)

        self.cumulative = cumulative

    def sample(self):
        """
        Perform the sampling

        Returns
        -------
        The sampled item.
        """
        r = self.rng.random()
        i = bisect.bisect_left(self.cumulative, r)
        return self.items[i]
