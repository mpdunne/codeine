from typing import Dict, Optional

from codeine.sequence.graph import CodonGraph, CodonNode


class SequenceSpace:
    """
    Class representing coding sequence space, for sampling and mutating CDS coding sequences.
    """

    def __init__(
            self,
            aa_seq: str,
            codon_restrictions: Optional[Dict[int, str]] = None,
            flank_l: str = '',
            flank_r: str = '',
    ) -> None:
        """
        Constructor for the SequenceSpace class.

        Parameters
        ----------
        aa_seq
            The amino acid sequence.
        codon_restrictions
            Any codon restrictions in the format e.g. {4: 'TCC'} or {5: ['AGT', 'AGC']}
        flank_l
            The context sequence to the left of the coding sequence
        flank_r
            The context sequence to the right of the coding sequence
        """
        self.graph = CodonGraph(
            aa_seq,
            codon_restrictions=codon_restrictions,
            flank_l=flank_l,
            flank_r=flank_r,
        )

    def sample(self, include_context: bool = False) -> str:
        node = self.graph.initial_node
        sequence = []

        while node.transitions:
            if isinstance(node, CodonNode):
                emitted = node.sample_codon()
                sequence.append(emitted)
                node = node.transitions[emitted]
            else:
                if include_context:
                    sequence.append(node.sequence)
                node = next(iter(node.transitions.values()))

        if include_context:
            sequence.append(node.sequence)

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
