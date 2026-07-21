from typing import Sequence


class Node:
    """
    Basic CodonGraph node.
    """

    def __init__(self, pos: int = None, id: str = None) -> None:
        """
        Constructor for the Node class.

        Parameters
        ----------
        pos
            The graph position of this Node.
        id
            The ID of this node.

        """
        self.pos = pos
        self.id = id

        self.parents = set()
        self.transitions = {}


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
        super().__init__(pos=pos, id=f'context-{pos}')

        # Basic info.
        self.sequence = sequence

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
        super().__init__(pos=pos, id=f'{aa}{pos}')

        # Basic info. Positioning is 1-based.
        self.aa = aa

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

    def __init__(self, pos: int) -> None:
        """
        Constructor for the EndNode class.

        pos
            The position in the graph. This should be len(aa_seq) + 2.
        """
        super().__init__(pos=pos, id='end')

    def __repr__(self) -> str:
        return (f'EndNode(id={self.id})')
