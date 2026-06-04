import uuid

from codeine.translation.tables import CodonTable
from codeine.utils.sampling import Sampler

from typing import Dict, List, Optional, Sequence, Union


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
        self.id = f"context-{uuid.uuid4().hex[:8]}"


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
        self.id = f"{aa}{pos}-{uuid.uuid4().hex[:8]}"

        # Initialise the basic attributes.
        self.codons = codons
        self.probabilities = [1] * len(self.codons)

        # Create the sampler.
        self.sampler = Sampler(self.codons, self.probabilities)

        # Option to pin (temporarily fix) a specific codon.
        self.pinned_codon = None

    def pin_codon(self, codon: str) -> None:
        """
        Pin (temporarily fix) a codon for this node.

        Parameters
        ----------
        codon
            The codon to pin.
        """
        codon = codon.upper()
        if codon not in self.codons:
            raise ValueError(f'Pinned codon {codon} is not a valid codon for position {self.pos}')
        self.pinned_codon = codon

    def unpin_codon(self) -> None:
        """
        Unpin any codons that are pinned on this node.
        """
        self.pinned_codon = None

    def sample_codon(self) -> str:
        """
        Sample a codon according to the stored weights.

        Returns
        -------
        A sampled codon.
        """
        if self.pinned_codon is not None:
            return self.pinned_codon
        return self.sampler.sample()


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
        self.id = f"end-{uuid.uuid4().hex[:8]}"


class CodonGraph:
    """
    Class representing a graph of codon nodes.
    """

    def __init__(
        self,
        aa_seq: str,
        codon_restrictions: Optional[Dict[int, CodonRestriction]] = None,
        context_l: str = '',
        context_r: str = '',
    ) -> None:
        if len(aa_seq) == 0:
            raise ValueError('Please provide non-empty sequence!')

        self.aa_seq = aa_seq.upper()
        self.codon_restrictions = codon_restrictions or {}
        self.ct = CodonTable()

        self.context_l = context_l.upper()
        self.context_r = context_r.upper()

        self.left_context_node = None
        self.right_context_node = None
        self.end_node = None

        self.codon_nodes = []
        self.codon_nodes_by_pos = {pos: [] for pos in range(1, len(self.aa_seq) + 1)}

        self.initial_node = None
        self.final_node = None

        self.node_counts = {}
        self.descendants_by_node = {}
        self.n_valid_sequences = 0

        self._validate_codon_restrictions()
        self.initialise_graph()

    def _validate_codon_restrictions(self) -> None:
        """
        Check the inputted codon restrictions make sense!
        """
        for pos, codon_restriction in self.codon_restrictions.items():
            if pos < 1 or pos > len(self.aa_seq):
                raise ValueError(f"Codon restriction position {pos} is out of range.")

            if isinstance(codon_restriction, str):
                codons = [codon_restriction]
            else:
                codons = list(codon_restriction)

            if len(codons) == 0:
                raise ValueError(f"Codon restriction at position {pos} cannot be empty.")

            codons = [codon.upper() for codon in codons]

            aa = self.aa_seq[pos - 1]
            allowed_codons = self.ct.aa_to_codons[aa]

            for codon in codons:
                if codon not in allowed_codons:
                    raise ValueError(f'Codon {codon} is not valid for amino acid {aa} at position {pos}.')

            self.codon_restrictions[pos] = codons

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
                codons = self.ct.aa_to_codons[aa]

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

        self.compile()

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

    def pin_codons(self, pinned_codons: Dict[int, str]) -> None:
        """
        Pin (temporarily fix) a codon in the codon graph.

        Parameters
        ----------
        pinned_codons
            A dict specifying which codons to pin, by pos: codon.
        """
        for pos, codon in pinned_codons.items():
            if pos not in self.codon_nodes_by_pos:
                raise ValueError(f'Pinned codon position {pos} is out of range.')

            for node in self.codon_nodes_by_pos[pos]:
                node.pin_codon(codon)

        self.compile()

    def unpin_codons(self, positions: List[int]) -> None:
        """
        Unpin codon nodes by pos.

        Parameters
        ----------
        positions
            A list of positions.
        """
        for pos in positions:
            if pos not in self.codon_nodes_by_pos:
                raise ValueError(f'Pinned codon position {pos} is out of range.')

            for node in self.codon_nodes_by_pos[pos]:
                node.unpin_codon()

        self.compile()

    def clear_pins(self) -> None:
        """
        Remove all codon pins from this graph.
        """
        for node in self.codon_nodes:
            node.unpin_codon()

        self.compile()

    @property
    def pinned_codons(self) -> Dict[int, str]:
        """
        Return a list of all codons that are pinned in this graph, and their codon values.

        Returns
        -------
        Dict keyed by pos and with codon values.
        """
        pinned = {}

        for pos, nodes in self.codon_nodes_by_pos.items():
            pinned_node_codons = {
                node.pinned_codon
                for node in nodes
                if node.pinned_codon is not None
            }

            if len(pinned_node_codons) == 1:
                pinned[pos] = next(iter(pinned_node_codons))
            elif len(pinned_node_codons) > 1:
                raise RuntimeError(
                    f'Position {pos} has inconsistent pinning across nodes.'
                )

        return pinned

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

    def _update_descendant_counts(self) -> None:
        """
        Calculate descendant counts for each node and outgoing transition.
        Assumes a strictly layered graph:

            left context -> codon positions -> right context -> end

        where all transitions move exactly one layer forwards.
        """
        node_counts = {self.end_node: 1}
        descendants_by_node = {}

        # Right context -> end
        transition_counts = {
            emitted: node_counts[child]
            for emitted, child in self.right_context_node.transitions.items()
        }
        descendants_by_node[self.right_context_node] = transition_counts
        node_counts[self.right_context_node] = sum(transition_counts.values())

        # Codon positions, backwards
        for pos in range(len(self.aa_seq), 0, -1):
            for node in self.codon_nodes_by_pos[pos]:
                transition_counts = {
                    emitted: node_counts[child]
                    for emitted, child in node.transitions.items()
                }

                descendants_by_node[node] = transition_counts
                node_counts[node] = sum(transition_counts.values())

        # Left context -> first codon position
        transition_counts = {
            emitted: node_counts[child]
            for emitted, child in self.left_context_node.transitions.items()
        }
        descendants_by_node[self.left_context_node] = transition_counts
        node_counts[self.left_context_node] = sum(transition_counts.values())

        self.node_counts = node_counts
        self.descendants_by_node = descendants_by_node
        self.n_valid_sequences = node_counts[self.initial_node]

    def compile(self) -> None:
        """
        Calculate all graph properties that are derived from its structure.
        Remember to do this after editing the graph!
        """
        # Calculate descendant counts!
        self._update_descendant_counts()
