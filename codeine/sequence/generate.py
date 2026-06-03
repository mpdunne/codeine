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
        """
        Pin (temporarily fix) a codon in the codon graph.

        Parameters
        ----------
        pinned_codons:
            A dict specifying which codons to pin, by pos: codon
        """
        self.graph.pin_codons(pinned_codons)

    def unpin_codons(self, positions):
        """
        Unpin codon nodes by pos.

        Parameters
        ----------
        positions:
            A list of positions
        """
        self.graph.unpin_codons(positions)

    def clear_pins(self):
        """
        Remove all codon pins from the generator.
        """
        self.graph.clear_pins()
