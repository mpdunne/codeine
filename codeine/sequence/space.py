from typing import Dict, Optional, Sequence

from codeine.sequence.graph import CodonGraph, CodonNode
from codeine.translation.tables import CodonTable


class SequenceSpace:
    """
    Class representing coding sequence space, for sampling and mutating CDS coding sequences.
    """

    def __init__(
            self,
            aa_seq: str,
            codon_restrictions: Optional[Dict[int, str]] = None,
            codon_table: CodonTable = None,
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
        codon_table
            The codon table to use. Leave blank to use standard table.
        context_l
            The context sequence to the left of the coding sequence
        context_r
            The context sequence to the right of the coding sequence
        """
        self.graph = CodonGraph(
            aa_seq,
            codon_restrictions=codon_restrictions,
            codon_table=codon_table,
            context_l=context_l,
            context_r=context_r,
        )

        self.user_pins = {}
        self.mutation_pins = {}

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

    def _update_graph_pins(self) -> None:
        """
        Update all active pins on the underlying graph.
        """
        self.graph.clear_pins()

        self.graph.pin_codons({
            **self.user_pins,
            **self.mutation_pins,
        })

    def pin_codons(self, pinned_codons):
        """
        Pin (temporarily fix) a codon in the codon graph.

        Parameters
        ----------
        pinned_codons:
            A dict specifying which codons to pin, by pos: codon
        """
        self.user_pins.update(pinned_codons)
        self._update_graph_pins()

    def unpin_codons(self, positions):
        """
        Unpin codon nodes by pos.

        Parameters
        ----------
        positions:
            A list of positions
        """
        for pos in positions:
            self.user_pins.pop(pos, None)

        self._update_graph_pins()

    def clear_pins(self):
        """
        Remove all codon pins from the generator.
        """
        self.user_pins.clear()
        self._update_graph_pins()

    def contains(self, seq: str) -> bool:
        """
        Check whether a DNA sequence is contained in this sequence space.

        Parameters
        ----------
        seq
            The sequence to check

        Returns
        -------
        True if and only if the sequence is contained in this sequence space.
        """
        return self.graph.contains(seq)

    def enter_mutation_mode(self, seq: str, positions: Sequence[int]) -> None:
        """
        Enter mutation mode, fixing the sequence on all but the specified positions

        Parameters
        ----------
        seq
            The sequence to mutate.
        positions
            The positions that are allowed to vary.
        """
        seq = seq.upper()

        if not self.contains(seq):
            raise ValueError("Parent sequence is not contained in this sequence space.")

        positions = set(positions)
        mutation_pins = {}

        for pos in range(1, len(self.graph.aa_seq) + 1):
            if pos not in positions:
                start = (pos - 1) * 3
                mutation_pins[pos] = seq[start:start + 3]

        self.mutation_pins = mutation_pins
        self._update_graph_pins()

    def exit_mutation_mode(self) -> None:
        """
        Clear all mutation restrictions and sample normally.
        """
        self.mutation_pins.clear()
        self._update_graph_pins()
