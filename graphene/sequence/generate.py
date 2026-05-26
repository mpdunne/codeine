import random

from graphene.sequence.graph import CodonGraph


class SequenceGenerator:
    """
    Basic coding sequence generator, using a coding sequence graph.
    """

    def __init__(self, aa_seq: str):
        """
        Constructor for the SequenceGenerator class.

        Parameters
        ----------
        aa_seq:
            The amino acid sequence.
        """
        self.graph = CodonGraph(aa_seq)

    def generate(self):
        initial_node = self.graph.initial_node
        node = initial_node
        sequence = ''
        while True:
            codon = random.choice(node.codons)
            sequence += codon
            if node.terminal:
                return sequence
            else:
                next_node = node.transitions[codon]
                node = next_node