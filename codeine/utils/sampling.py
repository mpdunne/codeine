import bisect
import random

from typing import Any, Sequence, Union


Seedable = Union[None, int, float, str, bytes, bytearray]


class Sampler:
    """
    A precomputed sampler to speed up weighted sampling.
    """
    def __init__(self, items: Sequence[Any], weights: Sequence[Union[int, float]] = None, seed: Seedable = None):
        """
        Constructor for the Sampler class.

        Parameters
        ----------
        items
            The items from which to sample.
        weights
            The weights assigned to the items.
        seed
            The random seed to use.
        """
        if len(items) == 0:
            raise ValueError('Items cannot be empty.')

        if weights is None:
            weights = [1] * len(items)

        if len(items) != len(weights):
            raise ValueError('Items and weights must have same length.')

        if any(weight < 0 for weight in weights):
            raise ValueError('Weights cannot be negative.')

        total = sum(weights)
        if total <= 0:
            raise ValueError('Weights must sum to a positive number')

        self.items = tuple(items)
        self._single = len(items) == 1

        if not self._single:
            self._rng = random.Random(seed)

            cumulative = []
            running = 0

            for weight in weights:
                running += weight
                cumulative.append(running / total)

            self._cumulative = cumulative

    def sample(self):
        """
        Sample the items according to the stored weights.
        If there is only one item, just return that.

        Returns
        -------
        The sampled item.
        """
        if self._single:
            return self.items[0]
        else:
            r = self._rng.random()
            i = bisect.bisect_left(self._cumulative, r)
            return self.items[i]
