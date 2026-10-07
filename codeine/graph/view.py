import random

from typing import Dict, Iterator, List, Optional, Sequence, Tuple, Union

from codeine.constraints.base import Constraint
from codeine.graph.base import CodonGraph, CodonRestriction
from codeine.graph.compiler import ViewCompiler
from codeine.graph.compiled import CompiledView
from codeine.graph.factorised import FactorisedCompiler
from codeine.translation.tables import TranslationTable
from codeine.translation.weights import CodonWeights
from codeine.utils.sampling import Seedable
from codeine.utils.tuples import tuplify


# Sentinels for compile status, indicating what, if
# anything, needs compiling/recompiling.
COMPILED = 0
COMPILE_SHALLOW = 1
COMPILE_EXTEND = 2
COMPILE_DEEP = 3


class CodonGraphView:
    """
    View of a codon graph with optional constraints and temporary codon pins.

    Holds configuration and manages lazy compilation. Sequence queries and sampling
    are delegated to the compiled result.
    """

    def __init__(self,
                 graph: CodonGraph,
                 *,
                 constraints: Optional[Union[Constraint, Sequence[Constraint]]] = None,
                 weights: Optional[CodonWeights] = None,
                 seed: Seedable = None,
                 compiler: str = 'flat',
                 ) -> None:
        """
        Constructor for the CodonGraphView

        Parameters
        ----------
        graph
            The underlying codon graph.
        constraints
            Any constraint trackers that we wish to use when traversing coding space.
        weights
            The codon weights to use when sampling.
        seed
            Seed used to initialise a random number generator.
        compiler
            Compilation engine, fixed at construction: ``'flat'`` (default) or
            ``'factorised'``.
        """

        if compiler not in ('flat', 'factorised'):
            raise ValueError(f"Unknown compiler {compiler!r}; supported compilers: 'flat', 'factorised'.")
        self._compiler = compiler

        self.graph = graph
        self.pinned_codons: Dict[int, CodonRestriction] = {}
        self.constraints: Tuple[Constraint, ...] = tuplify(constraints, Constraint)

        if weights is None:
            weights = CodonWeights.uniform(table=self.graph.tt)
        else:
            weights = weights.for_table(self.graph.tt)

        self._codon_weights = weights
        self._rng = random.Random(seed)

        self._compiled: Optional[CompiledView] = None
        self._compile_status = COMPILE_DEEP
        self._pending_constraints: Tuple[Constraint, ...] = ()

    def __getitem__(self, index: Union[int, slice]) -> Union[str, List[str]]:
        """
        Return one valid sequence, or a list of valid sequences for a slice.

        Parameters
        ----------
        index
            Zero-based sequence index, or slice of sequence indices.

        Returns
        -------
        str or list of str
            The indexed valid coding sequence, or a list of valid coding sequences.
        """
        if isinstance(index, slice):
            return self.sequences_at(index)

        elif isinstance(index, int):
            return self.sequence_at(index)

        else:
            raise ValueError(f'Invalid index: {index}')

    def __iter__(self) -> Iterator[str]:
        """
        Iterate over all valid sequences in this graph view.

        Yields
        ----------
        All valid sequences in the graph view, in order.
        """
        yield from self.enumerate()

    def __contains__(self, seq: str) -> bool:
        """
        Does the given seq exist in this space?

        Returns
        ----------
        True if and only if this is a valid sequence in this space.
        """
        return self.contains(seq)

    def contains(self, seq: str) -> bool:
        """
        Check whether a coding sequence is contained in this view.

        Parameters
        ----------
        seq
            The sequence to check

        Returns
        -------
        True if and only if the sequence is contained in this coding space.
        """
        if self._compile_status:
            self.compile()

        return self._compiled.contains(seq)

    def sample(self, n: Optional[int] = None) -> Union[str, List[str]]:
        """
        Sample one or more coding sequences from this graph view.

        Parameters
        ----------
        n
            Number of sequences to sample. If omitted, return a single sequence.

        Returns
        -------
        str or list of str
            One sampled coding sequence, or a list of sampled coding sequences.
        """
        if self._compile_status:
            self.compile()

        return self._compiled.sample(self._rng, n=n)

    def enumerate(self) -> Iterator[str]:
        """
        Enumerate all valid sequences in this view.

        Yields
        ------
        str
            All valid coding sequences, one by one.
        """
        if self._compile_status:
            self.compile()

        yield from self._compiled.enumerate()

    def enumerate_range(self, start: int = 0, stop: Optional[int] = None) -> Iterator[str]:
        """
        Enumerate valid sequences from start up to, but not including, stop.

        Parameters
        ----------
        start
            The zero-based start from which to begin enumeration
        stop
            The zero-based enumeration stop.

        Yields
        -------
        str
            Sequences in the range, one by one.
        """
        if self._compile_status:
            self.compile()

        yield from self._compiled.enumerate_range(start, stop)

    def sequence_at(self, index: int) -> str:
        """
        Return the valid sequence at a given index.

        Parameters
        ----------
        index
            Zero-based sequence index.

        Returns
        -------
        str
            The indexed valid coding sequence.
        """
        if self._compile_status:
            self.compile()

        return self._compiled.sequence_at(index)

    def sequences_at(self, index_slice: slice) -> List[str]:
        """
        Return valid sequences from a slice.
        """
        if self._compile_status:
            self.compile()

        return self._compiled.sequences_at(index_slice)

    def copy(self) -> 'CodonGraphView':
        """
        Copy this view and all its constraints and attributes.

        Returns
        -------
        A copy of the view.
        """
        view = self.graph.view(compiler=self.compiler)
        view._rng.setstate(self._rng.getstate())
        view._codon_weights = self._codon_weights

        view.pinned_codons = self.pinned_codons.copy()
        view.constraints = self.constraints

        view._compiled = self._compiled
        view._compile_status = self._compile_status
        view._pending_constraints = self._pending_constraints

        return view

    def compile(self) -> None:
        """
        Compile or recompile the view as required.
        """
        if self._compile_status == COMPILED:
            return

        compiler_type = {'flat': ViewCompiler, 'factorised': FactorisedCompiler}[self.compiler]
        compiler = compiler_type(self)

        if self._compiled is None or self._compile_status == COMPILE_DEEP:
            compiled = compiler.compile()

        elif self._compile_status == COMPILE_EXTEND:
            compiled = compiler.extend(self._compiled, self._pending_constraints)

        else:
            compiled = compiler.compile_shallow(self._compiled)

        self._compiled = compiled
        self._compile_status = COMPILED
        self._pending_constraints = ()

    def _update_compile_status(self, status: int) -> None:
        """
        Mark this view as requiring at least the specified compile phase.
        """
        self._compile_status = max(self._compile_status, status)

    def pin_codons(self, pinned_codons: Dict[int, CodonRestriction]) -> None:
        """
        Pin (temporarily fix) a codon in this codon graph view

        Parameters
        ----------
        pinned_codons
            A dict specifying which codons to pin, by pos: codon.
        """
        pinned_codons = self.graph.validate_fixed_codons(pinned_codons)
        self.pinned_codons.update(pinned_codons)
        self._update_compile_status(COMPILE_SHALLOW)

    def unpin_codons(self, positions: Union[int, Sequence[int]]) -> None:
        """
        Unpin codon nodes by pos.

        Parameters
        ----------
        positions
            A list of positions to unpin.
        """
        for pos in tuplify(positions, int):
            if pos < 1 or pos > len(self.graph.codon_nodes):
                raise ValueError(f'Pinned codon position {pos} is out of range.')

            self.pinned_codons.pop(pos, None)

        self._update_compile_status(COMPILE_SHALLOW)

    def set_pinned_codons(self, pinned_codons: Dict[int, CodonRestriction]) -> None:
        """
        Pin (temporarily fix) a specified group codons, leaving all others unpinned.

        Parameters
        ----------
        pinned_codons:
            A dict specifying which codons to pin, by pos: codon
        """
        pinned_codons = self.graph.validate_fixed_codons(pinned_codons)
        self.pinned_codons = dict(pinned_codons)
        self._update_compile_status(COMPILE_SHALLOW)

    def clear_pins(self) -> None:
        """
        Remove all codon pins from this graph view
        """
        self.pinned_codons.clear()
        self._update_compile_status(COMPILE_SHALLOW)

    def add_constraints(self, constraints: Union[Constraint, Sequence[Constraint]]) -> None:
        """
        Add one or more constraints to this view.

        Parameters
        ----------
        constraints
            Constraint or constraints to add.
        """
        constraints = tuplify(constraints, Constraint)

        if not constraints:
            return

        self.constraints += constraints
        self._pending_constraints += constraints
        self._update_compile_status(COMPILE_EXTEND)

    def set_constraints(self, constraints: Union[Constraint, Sequence[Constraint]]) -> None:
        """
        Set the constraints for this view.

        Parameters
        ----------
        constraints
            Constraints to apply during graph traversal.
        """
        self.constraints = tuplify(constraints, Constraint)
        self._pending_constraints = ()
        self._update_compile_status(COMPILE_DEEP)

    def clear_constraints(self) -> None:
        """
        Remove all constraints from this view.
        """
        self.set_constraints(())

    def set_weights(self, weights: Optional[CodonWeights] = None) -> None:
        """
        Set the codon weights used for sampling.

        Parameters
        ----------
        weights
            Codon weights to use. If omitted, use uniform weights.
        """
        if weights is None:
            weights = CodonWeights.uniform(table=self.graph.tt)
        else:
            weights = weights.for_table(self.graph.tt)

        self._codon_weights = weights
        self._update_compile_status(COMPILE_SHALLOW)

    def clear_weights(self) -> None:
        """
        Reset sampling to uniform codon weights.
        """
        self.set_weights()

    @property
    def compiler(self) -> str:
        """
        The compilation engine selected at construction.
        """
        return self._compiler

    @property
    def aa_seq(self) -> str:
        """
        The amino acid sequence.

        Returns
        -------
        The aa seq.
        """
        return self.graph.aa_seq

    @property
    def translation_table(self) -> TranslationTable:
        """
        The translation table used by the codon graph.
        """
        return self.graph.tt

    @property
    def codon_weights(self) -> CodonWeights:
        """
        The codon weights used by the codon graph.
        """
        return self._codon_weights

    @property
    def fixed_codons(self) -> Dict[int, CodonRestriction]:
        """
        Any hard-fixed codon restrictions on the codon graph.
        """
        return self.graph.fixed_codons

    @property
    def context_l(self) -> str:
        """
        The left context sequence.
        """
        return self.graph.context_l

    @property
    def context_r(self) -> str:
        """
        The right context sequence.
        """
        return self.graph.context_r

    @property
    def n_valid_sequences(self) -> int:
        """
        Number of valid coding sequences in this view given all constraints.
        """
        if self._compile_status:
            self.compile()

        return self._compiled.n_valid_sequences
