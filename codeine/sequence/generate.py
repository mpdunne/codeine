from typing import Dict, Optional

from codeine.sequence.graph import CodonGraph


class SequenceGenerator:
    """
    Basic coding sequence generator, using a coding sequence graph.
    """

    def __init__(self,
                 aa_seq: str,
                 codon_restrictions: Optional[Dict[int, str]] = None,
                 ) -> None:
        """
        Constructor for the SequenceGenerator class.

        Parameters
        ----------
        aa_seq:
            The amino acid sequence.
        """
        self.graph = CodonGraph(aa_seq, codon_restrictions=codon_restrictions)

    def generate(self) -> str:
        """
        Generate a DNA coding sequence for the given amino acid sequence!

        Returns
        -------
        A valid coding sequence.
        """
        initial_node = self.graph.initial_node
        node = initial_node
        codons = []
        while True:
            codon = node.sample_codon()
            codons.append(codon)
            if node.terminal:
                return ''.join(codons)
            else:
                next_node = node.transitions[codon]
                node = next_node

    def pin_codons(self, pinned_codons):
        self.graph.pin_codons(pinned_codons)

    def unpin_codons(self, positions):
        self.graph.unpin_codons(positions)

    def clear_pins(self):
        self.graph.clear_pins()
