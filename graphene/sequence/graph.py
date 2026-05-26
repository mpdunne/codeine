from graphene.translation.tables import CodonTable


class CodonNode:
    def __init__(self, pos, aa):
        self.pos = pos
        self.aa = aa
        self.codons = {}


class CodonGraph:
    def __init__(self, aa_seq):
        self.aa_seq = aa_seq.upper()
        self.ct = CodonTable()

        self.nodes = None
        self.initialise_graph()

    def initialise_graph(self):
        nodes = []

        for pos, aa in enumerate(self.aa_seq):
            node = CodonNode(pos, aa)

            codons = self.ct.aa_to_codons[aa]
            node.codons = {codon: None for codon in codons}
            nodes.append(node)

        for i in range(1, len(nodes)):
            previous = nodes[i - 1]
            current = nodes[i]
            previous.codons = {codon: current for codon in previous.codons}

        self.nodes = nodes
