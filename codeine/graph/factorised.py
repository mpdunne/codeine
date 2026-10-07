import random

from dataclasses import dataclass
from typing import Iterator, List, Optional, Union

from codeine.graph.base import CodonGraph
from codeine.graph.factor_count import ComponentModelCounter
from codeine.graph.factor_sampling import FactorSampler


@dataclass(frozen=True)
class FactorisedCompiledView:
    """
    A compiled factor model with the same sequence operations as a flat view.

    The counter owns validity and ordering; the sampler owns weights and its
    reusable plan. Neither stores the view's random generator or mutable pins.
    """
    graph: CodonGraph
    counter: ComponentModelCounter
    sampler: FactorSampler

    @property
    def n_valid_sequences(self) -> int:
        return self.counter.count()

    def contains(self, seq: str) -> bool:
        seq = self.graph.tt.normalise_sequence(seq)

        if len(seq) != len(self.graph.aa_seq) * 3:
            return False

        assignment = {pos: seq[(pos - 1) * 3:pos * 3] for pos in range(1, len(self.graph.aa_seq) + 1)}
        assignment[self.graph.left_context_node.pos] = self.graph.context_l
        assignment[self.graph.right_context_node.pos] = self.graph.context_r

        return self.counter.contains(assignment)

    def _sequence(self, assignment) -> str:
        return ''.join(assignment[pos] for pos in range(1, len(self.graph.aa_seq) + 1))

    def sample(self, rng: random.Random, n: Optional[int] = None) -> Union[str, List[str]]:
        if n is None:
            return self._sequence(self.sampler.sample(rng))

        return [self._sequence(assignment) for assignment in self.sampler.sample(rng, n)]

    def enumerate(self) -> Iterator[str]:
        return self.enumerate_range()

    def enumerate_range(self, start: int = 0, stop: Optional[int] = None) -> Iterator[str]:
        for assignment in self.counter.enumerate_range(start, stop):
            yield self._sequence(assignment)

    def sequence_at(self, index: int) -> str:
        return self._sequence(self.counter.assignment_at(index))

    def sequences_at(self, index_slice: slice) -> List[str]:
        start, stop, step = index_slice.indices(self.n_valid_sequences)

        if step == 1:
            return list(self.enumerate_range(start, max(start, stop)))

        return [self.sequence_at(index) for index in range(start, stop, step)]


class FactorisedCompiler:
    """
    Build factor models through the common constraint interface.

    Constraint types are interpreted by their own factors() methods. This
    compiler knows only domains, factors, pins and sampling weights.
    """

    def __init__(self, view) -> None:
        self.view = view
        self.graph = view.graph
        self.domains = {
            node.pos: tuple(
                choice for choice in node.transitions
                if node.pos not in view.pinned_codons or choice in view.pinned_codons[node.pos]
            )
            for node in self.graph.nodes
            if node is not self.graph.end_node
        }

    def _compiled(self, counter) -> FactorisedCompiledView:
        weights = {
            node.pos: {codon: self.view.codon_weights.weights[codon] for codon in node.codons}
            for node in self.graph.codon_nodes
        }
        counter.count()

        return FactorisedCompiledView(self.graph, counter, FactorSampler(counter, weights))

    def compile(self) -> FactorisedCompiledView:
        factors = []

        for constraint in self.view.constraints:
            constraint.link(self.graph)

            if constraint.is_trivial:
                continue

            local_factors = constraint.factors()

            if local_factors is None:
                raise ValueError(f'{type(constraint).__name__} does not support the factorised compiler.')

            factors.extend(local_factors)

        return self._compiled(ComponentModelCounter(self.domains, factors))

    def compile_shallow(self, compiled: FactorisedCompiledView) -> FactorisedCompiledView:
        """
        Reuse counts for weight changes; rebuild when pins change the domains.
        """
        if self.domains == compiled.counter.domains:
            return self._compiled(compiled.counter)

        return self.compile()

    def extend(self, compiled, constraints) -> FactorisedCompiledView:
        """
        Rebuild for added constraints, leaving the previous snapshot intact.
        """
        return self.compile()
