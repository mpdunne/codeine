import uuid

from graphene.translation.tables import CodonTable
from graphene.utils.sampling import Sampler

from typing import Dict, List, Optional


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

        # Initialise the sampler.
        self.codons = codons
        self.probabilities = [1] * len(codons)
        self.sampler = Sampler(self.codons, self.probabilities)

        # Graph stuff.
        self.transitions = {}
        self.parents = set()
        self.terminal = False

    def sample_codon(self) -> str:
        """
        Sample a codon according to the stored weights!

        Returns
        -------
        A sampled codon.
        """
        if self.sampler is None:
            raise ValueError(f"No sampler initialised for node {self.id}.")

        return self.sampler.sample()


class CodonGraph:
    """
    Class representing a graph of codon nodes.
    """

    def __init__(
        self,
        aa_seq: str,
        fixed_codons: Optional[Dict[int, str]] = None,
    ) -> None:
        if len(aa_seq) == 0:
            raise ValueError('Please provide non-empty sequence!')

        self.aa_seq = aa_seq.upper()
        self.fixed_codons = fixed_codons or {}
        self.ct = CodonTable()

        self.nodes = []
        self._initial_node = None

        self._validate_fixed_codons()
        self.initialise_graph()

    def _validate_fixed_codons(self) -> None:
        """
        Check the inputted fixed codons make sense!
        """
        for pos, codon in self.fixed_codons.items():
            if pos < 1 or pos > len(self.aa_seq):
                raise ValueError(f"Fixed codon position {pos} is out of range.")

            codon = codon.upper()
            aa = self.aa_seq[pos - 1]
            allowed_codons = self.ct.aa_to_codons[aa]

            if codon not in allowed_codons:
                raise ValueError(f'Codon {codon} is not valid for amino acid {aa} at position {pos}.')

            self.fixed_codons[pos] = codon

    def initialise_graph(self) -> None:
        """
        Initialise the codon graph.
        """
        nodes = []

        for ix, aa in enumerate(self.aa_seq):
            pos = ix + 1
            if pos in self.fixed_codons:
                codons = [self.fixed_codons[pos]]
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

    @property
    def initial_node(self) -> CodonNode:
        """
        The initial node in the graph. There's always only one initial node!

        Returns
        -------
        The initial node.
        """
        return self._initial_node
