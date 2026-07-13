import bisect
import random

from abc import ABC, abstractmethod
from typing import Any, Sequence, Optional, Union


Seedable = Union[int, float, str, bytes, bytearray]


class Sampler(ABC):
    """
    Generic abstract sampler class
    """
    @abstractmethod
    def sample(self):
        """
        Sample the item(s) according to whatever rules the inheriting
        class has decreed.

        Returns
        -------
        The sampled item.
        """
        pass


class SingletonSampler(Sampler):
    """
    A sample class that takes and returns a single item.
    """

    def __init__(
            self,
            item: Any,
    ):
        """
        Parameters
        ----------
        item
            The item to be "sampled".
        """
        super().__init__()
        self._item = item

    def sample(self):
        """
        "Sample" the item, i.e. return it.

        Returns
        -------
        The single item stored on this class.
        """
        return self._item


class RandomSampler(Sampler, ABC):
    """
    Base class for samplers that use an RNG to choose the returned items.
    """
    def __init__(
            self,
            items: Sequence[Any],
            seed: Optional[Seedable] = None,
            rng: Optional[random.Random] = None,
    ):
        """
        Parameters
        ----------
        items
            The items.
        seed
            Seed used to initialise a random number generator on this Sampler.
        rng
            Pre-constructed random number generator to use for sampling.
        """
        super().__init__()

        if len(items) == 0:
            raise ValueError('Items cannot be empty.')

        self._items = tuple(items)

        if seed is not None and rng is not None:
            raise ValueError('Provide either seed or rng, not both.')

        if rng is None:
            if seed is not None:
                rng = random.Random(seed)
            else:
                rng = random.Random()

        self._rng = rng


class UniformSampler(RandomSampler):
    """
    A uniform sampler.
    """
    def __init__(self,
                 items: Sequence[Any],
                 seed: Optional[Seedable] = None,
                 rng: Optional[random.Random] = None,
                 ):
        """
        Parameters
        ----------
        items
            The items from which to sample.
        seed
            Seed used to initialise a random number generator on this Sampler.
        rng
            Pre-constructed random number generator to use for sampling.
        """
        super().__init__(items=items, seed=seed, rng=rng)
        self._n = len(self._items)
        self._random = self._rng.random

    def sample(self):
        """
        Draw one of the stored items uniformly at random.

        Returns
        -------
        The sampled item.
        """
        # Note we're doing this because it's faster than random.sample(...), by about 1.5x.
        return self._items[int(self._random() * self._n)]


class WeightedSampler(RandomSampler):
    """
    A weighted sampler.
    """
    def __init__(self,
                 items: Sequence[Any],
                 weights: Sequence[Union[int, float]] = None,
                 seed: Optional[Seedable] = None,
                 rng: Optional[random.Random] = None,
                 ):
        """
        Parameters
        ----------
        items
            The items from which to sample.
        weights
            The weights assigned to the items.
        seed
            Seed used to initialise a random number generator on this Sampler.
        rng
            Pre-constructed random number generator to use for sampling.
        """
        super().__init__(items=items, seed=seed, rng=rng)

        if weights is None:
            weights = [1] * len(items)

        if len(items) != len(weights):
            raise ValueError('Items and weights must have same length.')

        if any(weight < 0 for weight in weights):
            raise ValueError('Weights cannot be negative.')

        total = sum(weights)
        if total <= 0:
            raise ValueError('Weights must sum to a positive number.')

        cumulative = []
        running = 0

        for weight in weights:
            running += weight
            cumulative.append(running / total)

        self._cumulative = cumulative

    def sample(self):
        """
        Sample the items according to the stored weights.

        Returns
        -------
        The sampled item.
        """
        r = self._rng.random()
        i = bisect.bisect_left(self._cumulative, r)
        return self._items[i]
