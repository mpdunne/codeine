from typing import Sequence


class Node:
    """
    Basic CodonGraph node.
    """

    def __init__(self) -> None:
        """
        Constructor for the Node class.
        """
        self.parents = set()
        self.transitions = {}
        self.pos = None


class ContextNode(Node):
    """
    Basic class representing a sequence context node on the codon graph.
    This refers to the sequence either to the left of or to the right of
    the coding sequence, and can be empty.
    """

    def __init__(self, pos: int, sequence: str) -> None:
        """
        Constructor for the ContextNode class.

        Parameters
        ----------
        pos
            The graph position. Left context is 0; right context is len(aa_seq) + 1.
        sequence
            The context sequence contained on this node.
        """
        super().__init__()

        # Basic info.
        self.pos = pos
        self.sequence = sequence

        # Set an ID for this node.
        self.id = f'context-{pos}'

    def __repr__(self) -> str:
        return (
            f'ContextNode('
            f'id={self.id}'
            f', pos={self.pos}'
            f')'
        )


class CodonNode(Node):
    """
    Basic class representing a codon node on the codon graph.
    """

    def __init__(self, pos: int, aa: str, codons: Sequence[str]) -> None:
        """
        Constructor for the CodonNode class.

        Parameters
        ----------
        pos
            The aa position. Positioning is 1-based.
        aa
            The aa identity.
        codons
            The possible codons for this node.
        """
        super().__init__()

        # Basic info. Positioning is 1-based.
        self.pos = pos
        self.aa = aa

        # Set an ID for this node.
        self.id = f'{aa}{pos}'

        # Initialise the basic attributes.
        self.codons = tuple(codons)

    def __repr__(self) -> str:
        codons = ','.join(self.codons)

        return (
            f'CodonNode('
            f'id={self.id}'
            f', pos={self.pos}'
            #f', codons=[{codons}]'
            f')'
        )


class EndNode(Node):
    """
    Final node for the codon graph.

    Marks successful completion of a graph walk.
    """

    def __init__(self) -> None:
        """
        Constructor for the EndNode class.
        """
        super().__init__()
        self.id = 'end'

    def __repr__(self) -> str:
        return (f'EndNode(id={self.id})')