import uuid

from graphene.translation.tables import CodonTable
from graphene.utils.sampling import Sampler

from typing import List


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
        # Basic info
        self.pos = pos
        self.aa = aa

        # Set an ID for this node.
        self.id = f"{aa}{pos}-{uuid.uuid4().hex[:8]}"

        # Initialise the sampler.
        self.codons = codons
        self.probabilities = [1] * len(codons)
        self.sampler = Sampler(self.codons, self.probabilities)

        # Graph stuff.
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
    def __init__(self, aa_seq: str) -> None:
        """
        Constructor for the CodonGraph class.

        Parameters
        ----------
        aa_seq:
            The amino acid sequence being represented.
        """
        if len(aa_seq) == 0:
            raise ValueError('Please provide non-empty sequence!')

        self.aa_seq = aa_seq.upper()
        self.ct = CodonTable()

        self.nodes = []
        self._initial_node = None
        self.initialise_graph()

    def initialise_graph(self) -> None:
        """
        Initialise the codon graph.

        Returns
        -------
        None
        """

        nodes = []
        for ix, aa in enumerate(self.aa_seq):
            pos = ix + 1
            codons = self.ct.aa_to_codons[aa]
            node = CodonNode(pos, aa, codons)
            nodes.append(node)

        self._initial_node = nodes[0]
        nodes[-1].terminal = True

        for i in range(1, len(nodes)):
            previous = nodes[i - 1]
            current = nodes[i]
            previous.transitions = {codon: current for codon in previous.codons}

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
