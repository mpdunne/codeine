from typing import Dict, Optional

from graphene.sequence.graph import CodonGraph


class SequenceGenerator:
    """
    Basic coding sequence generator, using a coding sequence graph.
    """

    def __init__(self,
                 aa_seq: str,
                 fixed_codons: Optional[Dict[int, str]] = None,
                 ) -> None:
        """
        Constructor for the SequenceGenerator class.

        Parameters
        ----------
        aa_seq:
            The amino acid sequence.
        """
        self.graph = CodonGraph(aa_seq, fixed_codons)

    def generate(self) -> str:
        """
        Generate a DNA coding sequence for the given amino acid sequence!

        Returns
        -------
        A valid coding sequence.
        """
        initial_node = self.graph.initial_node
        node = initial_node
        sequence = ''
        while True:
            codon = node.sample_codon()
            sequence += codon
            if node.terminal:
                return sequence
            else:
                next_node = node.transitions[codon]
                node = next_node