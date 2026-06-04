from typing import Dict, Optional

from codeine.sequence.graph import CodonGraph, CodonNode


class SequenceSpace:
    """
    Class representing coding sequence space, for sampling and mutating CDS coding sequences.
    """

    def __init__(self,
                 aa_seq: str,
                 codon_restrictions: Optional[Dict[int, str]] = None,
                 ) -> None:
        """
        Constructor for the SequenceSpace class.

        Parameters
        ----------
        aa_seq:
            The amino acid sequence.
        """
        self.graph = CodonGraph(aa_seq, codon_restrictions=codon_restrictions)

    def sample(self) -> str:
        """
        Sample a DNA coding sequence from this sequence space.

        Returns
        -------
        A valid coding sequence.
        """
        sequence = []
        node = self.graph.initial_node

        while node is not self.graph.final_node:

            if isinstance(node, CodonNode):
                emitted = node.sample_codon()
                sequence.append(emitted)
                node = node.transitions[emitted]
            else:
                node = next(iter(node.transitions.values()))

        return ''.join(sequence)

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
