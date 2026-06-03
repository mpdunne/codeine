import uuid

from codeine.translation.tables import CodonTable
from codeine.utils.sampling import Sampler

from typing import Dict, List, Optional, Sequence, Union


CodonRestriction = Union[str, Sequence[str]]


class CodonNode:
    """
    Basic class representing a node on the codon graph.
    """
    def __init__(self, pos: int, aa: str, codons: List[str]) -> None:
        """
        Constructor for the CodonNode class

        Parameters
        ----------
        pos
            The aa position (0-based).
        aa
            The aa identity.
        """
        # Basic info. Positioning is 1-based.
        self.pos = pos
        self.aa = aa

        # Set an ID for this node.
        self.id = f"{aa}{pos}-{uuid.uuid4().hex[:8]}"

        # Initialise the basic attributes
        self.codons = codons
        self.probabilities = [1] * len(codons)

        # Create the sampler
        self.sampler = Sampler(self.codons, self.probabilities)

        # Option to pin (temporarily fix) a specfic codon
        self.pinned_codon = None

        # Graph stuff.
        self.transitions = {}
        self.parents = set()
        self.terminal = False

    def pin_codon(self, codon : str):
        """
        Pin (temporarily fix) a codon for this node.

        Parameters
        ----------
        codon:
            The codon to pin
        """
        codon = codon.upper()
        if codon not in self.codons:
            raise ValueError(f'Pinned codon {codon} is not a valid codon for position {self.pos}')
        self.pinned_codon = codon

    def unpin_codon(self):
        """
        Unpin any codons that are pinned on this node.
        """
        self.pinned_codon = None

    def sample_codon(self) -> str:
        """
        Sample a codon according to the stored weights!

        Returns
        -------
        A sampled codon.
        """
        if self.pinned_codon is not None:
            return self.pinned_codon
        return self.sampler.sample()


class CodonGraph:
    """
    Class representing a graph of codon nodes.
    """

    def __init__(
        self,
        aa_seq: str,
        codon_restrictions: Optional[Dict[int, CodonRestriction]] = None,
    ) -> None:
        if len(aa_seq) == 0:
            raise ValueError('Please provide non-empty sequence!')

        self.aa_seq = aa_seq.upper()
        self.codon_restrictions = codon_restrictions or {}
        self.ct = CodonTable()

        self.nodes = []
        self.nodes_by_pos = {}
        self._initial_node = None

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
        nodes = []

        for ix, aa in enumerate(self.aa_seq):
            pos = ix + 1
            if pos in self.codon_restrictions:
                codons = self.codon_restrictions[pos]
            else:
                codons = self.ct.aa_to_codons[aa]

            node = CodonNode(pos, aa, codons)
            nodes.append(node)

        for i in range(1, len(nodes)):
            previous = nodes[i - 1]
            current = nodes[i]
            previous.transitions = {
                codon: current
                for codon in previous.codons
            }
            current.parents.add(previous)

        self._initial_node = nodes[0]
        nodes[-1].terminal = True
        self.nodes = nodes

        for node in nodes:
            self.nodes_by_pos.setdefault(node.pos, []).append(node)

    @property
    def initial_node(self) -> CodonNode:
        """
        The initial node in the graph. There's always only one initial node!

        Returns
        -------
        The initial node.
        """
        return self._initial_node

    def pin_codons(self, pinned_codons: Dict[int, str]) -> None:
        """
        Pin (temporarily fix) a codon in the codon graph.

        Parameters
        ----------
        pinned_codons:
            A dict specifying which codons to pin, by pos: codon
        """
        for pos, codon in pinned_codons.items():
            if pos not in self.nodes_by_pos:
                raise ValueError(f'Pinned codon position {pos} is out of range.')

            for node in self.nodes_by_pos[pos]:
                node.pin_codon(codon)

    def unpin_codons(self, positions: List[int]) -> None:
        """
        Unpin codon nodes by pos.

        Parameters
        ----------
        positions:
            A list of positions
        """
        for pos in positions:
            if pos not in self.nodes_by_pos:
                raise ValueError(f'Pinned codon position {pos} is out of range.')

            for node in self.nodes_by_pos[pos]:
                node.unpin_codon()

    def clear_pins(self) -> None:
        """
        Remove all codon pins from this graph.
        """
        for node in self.nodes:
            node.unpin_codon()

    @property
    def pinned_codons(self) -> Dict[int, str]:
        """
        Return a list of all codons that are pinned in this graph, and their codon values.

        Returns
        -------
        Dict keyed by pos and with codon values.
        """
        pinned = {}

        for pos, nodes in self.nodes_by_pos.items():
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
