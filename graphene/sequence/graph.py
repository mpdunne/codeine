import uuid

from graphene.translation.tables import CodonTable
from graphene.utils.sampling import Sampler


class CodonNode:
    """
    Basic class representing a node on the codon graph.
    """
    def __init__(self, pos: int, aa: str):
        """
        Constructor for the CodonNode class

        Parameters
        ----------
        pos
            The aa position (0-based).
        aa
            The aa identity.
        """
        self.pos = pos
        self.aa = aa

        # These will be set separately.
        self.codons = []
        self.probabilities = []
        self.transitions = {}
        self.terminal = False
        self.sampler = None

        # Set an ID for this node.
        self.id = f"{aa}{pos}-{uuid.uuid4().hex[:8]}"

    def sample_codon(self):
        if self.sampler is None:
            raise ValueError(f"No sampler initialised for node {self.id}.")

        return self.sampler.sample()


class CodonGraph:
    """
    Class representing a graph of codon nodes.
    """
    def __init__(self, aa_seq: str):
        if len(aa_seq) == 0:
            raise ValueError('Please provide non-empty sequence!')

        self.aa_seq = aa_seq.upper()
        self.ct = CodonTable()

        self.nodes = []
        self._initial_node = None
        self.initialise_graph()

    def initialise_graph(self):

        nodes = []
        for pos, aa in enumerate(self.aa_seq):
            node = CodonNode(pos, aa)
            node.codons = self.ct.aa_to_codons[aa]
            node.probabilities = [1] * len(node.codons)
            node.sampler = Sampler(node.codons, node.probabilities)
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
        return self._initial_node
