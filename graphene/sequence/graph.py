from graphene.translation.tables import CodonTable


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
        self.codons = []
        self.transitions = {}
        self.terminal = False


class CodonGraph:
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
