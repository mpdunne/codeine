from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from codeine.sequence.view import CodonGraphView

import uuid

from collections import Counter
from typing import Dict, List, Optional, Sequence, Union, Tuple

from codeine.sequence.display import format_banned_sequences, format_restrictions
from codeine.translation.tables import TranslationTable
from codeine.translation.weights import CodonWeights


CodonRestriction = Union[str, Sequence[str]]


class Node:
    """
    Basic CodonGraph node.
    """

    def __init__(self) -> None:
        """
        Constructor for the Node class.
        """
        self.parents = set()
        self.transitions = {}


class ContextNode(Node):
    """
    Basic class representing a sequence context node on the codon graph.
    This refers to the sequence either to the left of or to the right of
    the coding sequence, and can be empty.
    """

    def __init__(self, sequence: str) -> None:
        """
        Constructor for the ContextNode class.

        Parameters
        ----------
        sequence
            The context sequence emitted by this node.
        """
        super().__init__()

        # Basic info.
        self.sequence = sequence

        # Set an ID for this node.
        self.id = f'context-{uuid.uuid4().hex[:8]}'


class CodonNode(Node):
    """
    Basic class representing a codon node on the codon graph.
    """

    def __init__(self, pos: int, aa: str, codons: List[str]) -> None:
        """
        Constructor for the CodonNode class.

        Parameters
        ----------
        pos
            The aa position. Positioning is 1-based.
        aa
            The aa identity.
        codons
            The possible codons for this node.
        """
        super().__init__()

        # Basic info. Positioning is 1-based.
        self.pos = pos
        self.aa = aa

        # Set an ID for this node.
        self.id = f'{aa}{pos}-{uuid.uuid4().hex[:8]}'

        # Initialise the basic attributes.
        self.codons = codons


class EndNode(Node):
    """
    Final node for the codon graph.

    Marks successful completion of a graph walk.
    """

    def __init__(self) -> None:
        """
        Constructor for the EndNode class.
        """
        super().__init__()
        self.id = f'end-{uuid.uuid4().hex[:8]}'


class CodonGraph:
    """
    Class representing a graph of codon nodes.
    """

    def __init__(
        self,
        aa_seq: str,
        codon_restrictions: Optional[Dict[int, CodonRestriction]] = None,
        translation_table: TranslationTable = None,
        weights: CodonWeights = None,
        context_l: str = '',
        context_r: str = '',
    ) -> None:
        if len(aa_seq) == 0:
            raise ValueError('Please provide non-empty sequence!')

        self.aa_seq = aa_seq.upper()

        if translation_table is None:
            rna = weights.rna if weights is not None else False
            translation_table = TranslationTable(table_id=1, rna=rna)

        if weights is None:
            weights = CodonWeights.uniform(table=translation_table, rna=translation_table.rna)

        self.validate_codon_weights(weights, translation_table)

        self.tt = translation_table
        self.cw = weights

        self.codon_restrictions = {}
        self.codon_restrictions = self.validate_codon_restrictions(codon_restrictions)

        self.banned_sequences = []

        self.context_l = context_l.upper()
        self.context_r = context_r.upper()

        self.left_context_node = None
        self.right_context_node = None
        self.end_node = None

        self.codon_nodes = []
        self.codon_nodes_by_pos = {
            pos: []
            for pos in range(1, len(self.aa_seq) + 1)
        }

        self.initial_node = None
        self.final_node = None

        self.initialise_graph()

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

        if self.banned_sequences:
            lines += [
                'Banned sequences:',
                *format_banned_sequences(self.banned_sequences),
            ]

        return '\n'.join(lines)

    def validate_codon_restrictions(self, codon_restrictions: Dict[int, CodonRestriction]) -> Dict[int, List[str]]:
        """
        Check the inputted codon restrictions make sense!
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

            codons = [codon.upper() for codon in codons]

            aa = self.aa_seq[pos - 1]

            if pos in self.codon_restrictions:
                allowed_codons = [self.tt.normalise_codon(codon) for codon in self.codon_restrictions[pos]]
            else:
                allowed_codons = self.tt.aa_to_codons[aa]

            for codon in codons:
                if codon not in allowed_codons:
                    raise ValueError(f'Codon {codon} is not valid for amino acid {aa} at position {pos}.')

            normalised[pos] = codons

        return normalised

    @staticmethod
    def validate_codon_weights(
            weights: CodonWeights,
            translation_table: TranslationTable,
    ) -> None:
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

    def initialise_graph(self) -> None:
        """
        Initialise the codon graph.
        """
        left_context_node = ContextNode(self.context_l)
        right_context_node = ContextNode(self.context_r)
        end_node = EndNode()

        for ix, aa in enumerate(self.aa_seq):
            pos = ix + 1

            if pos in self.codon_restrictions:
                codons = self.codon_restrictions[pos]
            else:
                codons = self.tt.aa_to_codons[aa]

            node = CodonNode(pos, aa, codons)
            self.add_codon_node(node)

        # Left context -> first codon node
        left_context_node.transitions = {
            left_context_node.sequence: self.codon_nodes[0]
        }
        self.codon_nodes[0].parents.add(
            (left_context_node, left_context_node.sequence)
        )

        # Codon node -> next codon node
        for i in range(1, len(self.codon_nodes)):
            previous = self.codon_nodes[i - 1]
            current = self.codon_nodes[i]

            for codon in previous.codons:
                previous.transitions[codon] = current
                current.parents.add((previous, codon))

        # Last codon node -> right context
        last_codon_node = self.codon_nodes[-1]
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

        self.initial_node = left_context_node
        self.final_node = end_node

    def _find_matching_subpaths(self, sequence: str) -> List[Tuple[List[Tuple[Node, str]], int]]:
        """
        For a given sequence, find subpaths in the graph that match that sequence.
        Return each found subpath in the following format:

            (
                [
                    (node1, codon_1),
                    (node2, codon_2),
                    ...
                ]
                offset,  # Where the path starts relative to first node's codon choice.
            )

        Parameters
        ----------
        sequence
            The sequence to search for.

        Returns
        -------
        A tuple consisteing of a list of (node, codon) pairs, followed by the offset.
        """

        sequence = sequence.upper()

        if len(sequence) == 0:
            raise ValueError('Sequence cannot be empty.')

        matches = []
        candidate_matches = []

        # First, check which nodes we can start at.
        for node in self.nodes:
            if node is self.end_node:
                continue

            for choice, child in node.transitions.items():
                for offset in range(len(choice)):
                    choice_subsequence = choice[offset:]

                    if choice_subsequence.startswith(sequence):
                        # Bingo!
                        matches.append(([(node, choice)], offset))

                    elif sequence.startswith(choice_subsequence):
                        # Maybe bingo! Maygo!
                        candidate_matches.append((
                            [(node, choice)],
                            offset,
                            len(choice_subsequence),
                        ))

        def reinspect_candidate_matches(candidate_matches):
            reinspect = []

            for partial_path, offset, seen_length in candidate_matches:
                previous_node, previous_choice = partial_path[-1]
                node = previous_node.transitions[previous_choice]

                remaining_sequence = sequence[seen_length:]

                if remaining_sequence == '':
                    # Fantastic!
                    matches.append((partial_path, offset))
                    continue

                if node is self.final_node:
                    continue

                if isinstance(node, CodonNode):
                    choice_length = 3
                else:
                    choice_length = len(next(iter(node.transitions)))

                # Sneaky shortcut if we've crossed into the right context:
                if isinstance(node, CodonNode):
                    pos = node.pos

                    remaining_sequence_length = len(sequence) - seen_length
                    remaining_coding_length = 3 * (len(self.aa_seq) - pos + 1)

                    if remaining_sequence_length > remaining_coding_length:
                        sequence_end = sequence[seen_length + remaining_coding_length:]

                        if not self.context_r.startswith(sequence_end):
                            continue

                for choice, child in node.transitions.items():

                    if len(remaining_sequence) >= choice_length:

                        if remaining_sequence.startswith(choice):
                            # Keep going...
                            reinspect.append((partial_path + [(node, choice)], offset, seen_length + choice_length))

                        else:
                            # Hard luck this time.
                            continue

                    else:

                        if choice.startswith(remaining_sequence):
                            # Wahoo!
                            matches.append((partial_path + [(node, choice)], offset))

                        else:
                            # Hard luck this time.
                            continue

            return reinspect

        while candidate_matches:
            candidate_matches = reinspect_candidate_matches(candidate_matches)

        return matches

    @property
    def nodes(self) -> List[Node]:
        """
        All nodes in the graph, including context and end nodes.
        """
        return [
            self.left_context_node,
            *self.codon_nodes,
            self.right_context_node,
            self.end_node,
        ]

    def add_codon_node(self, node: CodonNode) -> None:
        """
        Add a codon node and update codon-node indexes.
        """
        self.codon_nodes.append(node)
        self.codon_nodes_by_pos.setdefault(node.pos, []).append(node)

    def remove_codon_node(self, node: CodonNode) -> None:
        """
        Remove a codon node and update codon-node indexes.
        """
        self.codon_nodes.remove(node)
        self.codon_nodes_by_pos[node.pos].remove(node)

    def view(self) -> 'CodonGraphView':
        """
        Return a constrained view over this graph.
        """
        from codeine.sequence.view import CodonGraphView
        return CodonGraphView(self)
