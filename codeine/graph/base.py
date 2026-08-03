from __future__ import annotations

from typing import TYPE_CHECKING, Dict, List, Optional, Sequence, Tuple, Union

from codeine.graph.nodes import CodonNode, ContextNode, EndNode, Node
from codeine.translation.tables import TranslationTable
from codeine.utils.display import format_restrictions
from codeine.utils.sampling import Seedable

if TYPE_CHECKING:
    from codeine.constraints.base import Constraint
    from codeine.graph.view import CodonGraphView
    from codeine.translation.weights import CodonWeights

CodonRestriction = Union[str, Sequence[str]]


class CodonGraph:
    """
    Class representing a graph of codon nodes.
    """

    def __init__(
        self,
        aa_seq: str,
        fixed_codons: Optional[Dict[int, CodonRestriction]] = None,
        translation_table: Optional[TranslationTable] = None,
        context_l: str = '',
        context_r: str = '',
    ) -> None:
        """
        Parameters
        ----------
        aa_seq
            The amino acid sequence.
        fixed_codons
            Any codon restrictions, for example fixed codons or subsets, in the form {3: 'AAA', 6: ['TTT', 'TTC']...}
        translation_table
            The translation table to use.
        context_l
            The left context sequence.
        context_r
            The right context sequence.
        """
        if len(aa_seq) == 0:
            raise ValueError('Please provide non-empty sequence!')

        if translation_table is None:
            translation_table = TranslationTable(table_id=1, rna=False)

        self.tt = translation_table

        self.aa_seq = aa_seq.upper()
        self.validate_aa_seq()

        self.fixed_codons = {}
        self.fixed_codons = self.validate_fixed_codons(fixed_codons)

        self.context_l = self.tt.normalise_sequence(context_l)
        self.context_r = self.tt.normalise_sequence(context_r)

        self.left_context_node = None
        self.right_context_node = None
        self.end_node = None

        self.codon_nodes: Tuple[CodonNode, ...] = ()

        self.initial_node = None
        self.final_node = None

        self._initialise_graph()

    def validate_aa_seq(self) -> None:
        """
        Check that all amino acids in the sequence are supported.
        """
        for pos, aa in enumerate(self.aa_seq, start=1):
            if aa not in self.tt.aa_to_codons:
                raise ValueError(f'Invalid amino acid {aa} at position {pos}.')

    def validate_fixed_codons(
            self,
            fixed_codons: Optional[Dict[int, CodonRestriction]],
    ) -> Dict[int, List[str]]:
        """
        Check the inputted restrictions make sense!
        """
        fixed_codons = fixed_codons or {}
        normalised = {}

        for pos, codon_restriction in fixed_codons.items():
            if pos < 1 or pos > len(self.aa_seq):
                raise ValueError(f'Restricted position {pos} is out of range.')

            if isinstance(codon_restriction, str):
                codons = [codon_restriction]
            else:
                codons = list(set(codon_restriction))

            if len(codons) == 0:
                raise ValueError(f'Codon restriction at position {pos} cannot be empty.')

            codons = [self.tt.normalise_sequence(codon) for codon in codons]

            aa = self.aa_seq[pos - 1]

            if pos in self.fixed_codons:
                allowed_codons = [self.tt.normalise_sequence(codon) for codon in self.fixed_codons[pos]]
            else:
                allowed_codons = self.tt.aa_to_codons[aa]

            for codon in codons:
                if codon not in allowed_codons:
                    raise ValueError(f'Codon {codon} is not valid for amino acid {aa} at position {pos}.')

            normalised[pos] = codons

        return normalised

    def codon_node_by_pos(self, pos: int) -> CodonNode:
        """
        Return the codon node at a given amino-acid position.

        Positioning is 1-based.
        """
        if pos < 1 or pos > len(self.codon_nodes):
            raise ValueError(f'Position {pos} is out of range.')

        return self.codon_nodes[pos - 1]

    def _initialise_graph(self) -> None:
        """
        Initialise the codon graph.
        """
        left_context_node = ContextNode(pos=0, sequence=self.context_l)
        right_context_node = ContextNode(pos=len(self.aa_seq) + 1, sequence=self.context_r)
        end_node = EndNode(pos=len(self.aa_seq) + 2)

        codon_nodes = []
        for ix, aa in enumerate(self.aa_seq):
            pos = ix + 1

            if pos in self.fixed_codons:
                codons = self.fixed_codons[pos]
            else:
                codons = self.tt.aa_to_codons[aa]

            node = CodonNode(pos=pos, aa=aa, codons=codons)
            codon_nodes.append(node)

        # Left context -> first codon node
        left_context_node.transitions = {
            left_context_node.sequence: codon_nodes[0]
        }
        codon_nodes[0].parents.add(
            (left_context_node, left_context_node.sequence)
        )

        # Codon node -> next codon node
        for i in range(1, len(codon_nodes)):
            previous = codon_nodes[i - 1]
            current = codon_nodes[i]

            for codon in previous.codons:
                previous.transitions[codon] = current
                current.parents.add((previous, codon))

        # Last codon node -> right context
        last_codon_node = codon_nodes[-1]
        for codon in last_codon_node.codons:
            last_codon_node.transitions[codon] = right_context_node
            right_context_node.parents.add((last_codon_node, codon))

        # Right context -> end
        right_context_node.transitions = {
            right_context_node.sequence: end_node
        }
        end_node.parents.add(
            (right_context_node, right_context_node.sequence)
        )

        self.left_context_node = left_context_node
        self.right_context_node = right_context_node
        self.end_node = end_node
        self.codon_nodes = tuple(codon_nodes)

        self.initial_node = left_context_node
        self.final_node = end_node

    @property
    def nodes(self) -> Tuple[Node, ...]:
        """
        All nodes in the graph, including context and end nodes.
        """
        return (
            self.left_context_node,
            *self.codon_nodes,
            self.right_context_node,
            self.end_node,
        )

    def view(self,
             *,
             constraints: Optional[Sequence[Constraint]] = None,
             weights: Optional[CodonWeights] = None,
             seed: Seedable = None,
             ) -> 'CodonGraphView':
        """
        Return

        Parameters
        ----------
        constraints
            Any constraint trackers that we wish to use when traversing coding space.
        weights
            The codon weights to use when sampling.
        seed
            Seed used to initialise a random number generator.

        Returns
        -------
        A constrained view over this graph.
        """
        from codeine.graph.view import CodonGraphView
        return CodonGraphView(self, seed=seed, weights=weights, constraints=constraints)
