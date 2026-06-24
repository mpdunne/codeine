from __future__ import annotations

import random

from collections import Counter
from typing import TYPE_CHECKING, Dict, List, Optional, Sequence, Tuple, Union

from codeine.graph.nodes import CodonNode, ContextNode, EndNode, Node
from codeine.translation.tables import TranslationTable
from codeine.translation.weights import CodonWeights
from codeine.utils.display import format_restrictions
from codeine.utils.sampling import Seedable

if TYPE_CHECKING:
    from codeine.graph.view import CodonGraphView

CodonRestriction = Union[str, Sequence[str]]


class CodonGraph:
    """
    Class representing a graph of codon nodes.
    """

    def __init__(
        self,
        aa_seq: str,
        codon_restrictions: Optional[Dict[int, CodonRestriction]] = None,
        translation_table: Optional[TranslationTable] = None,
        weights: Optional[CodonWeights] = None,
        context_l: str = '',
        context_r: str = '',
    ) -> None:
        if len(aa_seq) == 0:
            raise ValueError('Please provide non-empty sequence!')

        if translation_table is None:
            rna = weights.rna if weights is not None else False
            translation_table = TranslationTable(table_id=1, rna=rna)

        if weights is None:
            weights = CodonWeights.uniform(table=translation_table, rna=translation_table.rna)

        self.validate_codon_weights(weights, translation_table)

        self.tt = translation_table
        self.cw = weights

        self.aa_seq = aa_seq.upper()
        self.validate_aa_seq()

        self.codon_restrictions = {}
        self.codon_restrictions = self.validate_codon_restrictions(codon_restrictions)

        self.context_l = self.tt.normalise_sequence(context_l)
        self.context_r = self.tt.normalise_sequence(context_r)

        self.left_context_node = None
        self.right_context_node = None
        self.end_node = None

        self.codon_nodes: Tuple[CodonNode, ...] = ()

        self.initial_node = None
        self.final_node = None

        self._initialise_graph()

    def __repr__(self) -> str:
        molecule = 'RNA' if self.tt.rna else 'DNA'

        lines = [
            f'{type(self).__name__}',
            '',
            f'Translation table: {self.tt.table_id} ({self.tt.name})',
            f'Molecule type: {molecule}',
            '',
            f'Amino acid sequence ({len(self.aa_seq)} aa)',
            f'{self.aa_seq}',
            ''
        ]
        if self.codon_restrictions:
            lines += [
                'Codon restrictions:',
                *format_restrictions(
                    self.codon_restrictions,
                    label='restricted positions',
                ),
                '',
            ]

        return '\n'.join(lines)

    def validate_aa_seq(self) -> None:
        """
        Check that all amino acids in the sequence are supported.
        """
        for pos, aa in enumerate(self.aa_seq, start=1):
            if aa not in self.tt.aa_to_codons:
                raise ValueError(f'Invalid amino acid {aa} at position {pos}.')

    def validate_codon_restrictions(
            self,
            codon_restrictions: Optional[Dict[int, CodonRestriction]],
    ) -> Dict[int, List[str]]:
        """
        Check the inputted restrictions make sense!
        """
        codon_restrictions = codon_restrictions or {}
        normalised = {}

        for pos, codon_restriction in codon_restrictions.items():
            if pos < 1 or pos > len(self.aa_seq):
                raise ValueError(f'Restricted position {pos} is out of range.')

            if isinstance(codon_restriction, str):
                codons = [codon_restriction]
            else:
                codons = list(codon_restriction)

            if len(codons) == 0:
                raise ValueError(f'Codon restriction at position {pos} cannot be empty.')

            codons = [self.tt.normalise_sequence(codon) for codon in codons]

            aa = self.aa_seq[pos - 1]

            if pos in self.codon_restrictions:
                allowed_codons = [self.tt.normalise_sequence(codon) for codon in self.codon_restrictions[pos]]
            else:
                allowed_codons = self.tt.aa_to_codons[aa]

            for codon in codons:
                if codon not in allowed_codons:
                    raise ValueError(f'Codon {codon} is not valid for amino acid {aa} at position {pos}.')

            normalised[pos] = codons

        return normalised

    @staticmethod
    def validate_codon_weights(weights: CodonWeights, translation_table: TranslationTable) -> None:
        """
        Check that codon weights are compatible with the provided translation table.

        Parameters
        ----------
        weights
            The codon weights.
        translation_table
            The translation table.

        Raises
        -------
        Various errors if things aren't good.
        """
        if weights.rna != translation_table.rna:
            raise ValueError('Codon weights and translation table use different molecule types.')

        expected_codons = {
            aa: Counter(codons)
            for aa, codons in translation_table.aa_to_codons.items()
        }

        actual_codons = {
            aa: Counter(codons)
            for aa, codons in weights.aa_to_codons.items()
        }

        if actual_codons != expected_codons:
            raise ValueError('Codon weights and translation table do not match.')

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
        left_context_node = ContextNode(0, self.context_l)
        right_context_node = ContextNode(len(self.aa_seq) + 1, self.context_r)
        end_node = EndNode()

        codon_nodes = []
        for ix, aa in enumerate(self.aa_seq):
            pos = ix + 1

            if pos in self.codon_restrictions:
                codons = self.codon_restrictions[pos]
            else:
                codons = self.tt.aa_to_codons[aa]

            node = CodonNode(pos, aa, codons)
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
             seed: Optional[Seedable] = None,
             rng: Optional[random.Random] = None) -> 'CodonGraphView':
        """
        Return a constrained view over this graph.

        Parameters
        ----------
        seed
            Seed used to initialise the view's random number generator.
        rng
            Random number generator used by the view for sampling.
        """
        from codeine.graph.view import CodonGraphView
        return CodonGraphView(self, seed=seed, rng=rng)
