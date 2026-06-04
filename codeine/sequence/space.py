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
            context_l: str = '',
            context_r: str = '',
    ) -> None:
        """
        Constructor for the SequenceSpace class.

        Parameters
        ----------
        aa_seq
            The amino acid sequence.
        codon_restrictions
            Any codon restrictions in the format e.g. {4: 'TCC'} or {5: ['AGT', 'AGC']}
        context_l
            The context sequence to the left of the coding sequence
        context_r
            The context sequence to the right of the coding sequence
        """
        self.graph = CodonGraph(
            aa_seq,
            codon_restrictions=codon_restrictions,
            flank_l=context_l,
            flank_r=context_r,
        )

    def sample(self, include_context: bool = False) -> str:
        """
        Sample a DNA sequence from this sequence space.

        Parameters
        ----------
        include_context
            Whether to include the left and right context sequences.

        Returns
        -------
        A sampled sequence. By default, only the coding sequence is returned.
        """
        node = self.graph.initial_node
        sequence = []

        while node is not self.graph.final_node:
            if not node.transitions:
                raise RuntimeError(f"Reached non-final node {node.id} with no outgoing transitions.")

            if isinstance(node, CodonNode):
                emitted = node.sample_codon()
                sequence.append(emitted)

            else:
                emitted = node.sequence
                if include_context:
                    sequence.append(emitted)

            node = node.transitions[emitted]

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
