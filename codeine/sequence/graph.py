import uuid

from codeine.translation.tables import CodonTable
from codeine.utils.sampling import Sampler

from typing import Dict, Generator, List, Optional, Sequence, Union


CodonRestriction = Union[str, Sequence[str]]


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


class ContextNode(Node):
    """
    Basic class representing a sequence context node on the codon graph.
    This refers to the sequence either to the left of or to the right of
    the coding sequence, and can be empty.
    """

    def __init__(self, sequence: str) -> None:
        """
        Constructor for the ContextNode class.

        Parameters
        ----------
        sequence
            The context sequence emitted by this node.
        """
        super().__init__()

        # Basic info.
        self.sequence = sequence

        # Set an ID for this node.
        self.id = f"context-{uuid.uuid4().hex[:8]}"


class CodonNode(Node):
    """
    Basic class representing a codon node on the codon graph.
    """

    def __init__(self, pos: int, aa: str, codons: List[str]) -> None:
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
        self.id = f"{aa}{pos}-{uuid.uuid4().hex[:8]}"

        # Initialise the basic attributes.
        self.codons = codons
        self.probabilities = [1] * len(self.codons)


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
        self.id = f"end-{uuid.uuid4().hex[:8]}"


class CodonGraph:
    """
    Class representing a graph of codon nodes.
    """

    def __init__(
        self,
        aa_seq: str,
        codon_restrictions: Optional[Dict[int, CodonRestriction]] = None,
        codon_table: CodonTable = None,
        context_l: str = '',
        context_r: str = '',
    ) -> None:
        if len(aa_seq) == 0:
            raise ValueError('Please provide non-empty sequence!')

        self.aa_seq = aa_seq.upper()

        if codon_table is None:
            codon_table = CodonTable()
        self.ct = codon_table

        self.codon_restrictions = {}
        self.codon_restrictions = self.validate_codon_restrictions(codon_restrictions)

        self.context_l = context_l.upper()
        self.context_r = context_r.upper()

        self.left_context_node = None
        self.right_context_node = None
        self.end_node = None

        self.codon_nodes = []
        self.codon_nodes_by_pos = {pos: [] for pos in range(1, len(self.aa_seq) + 1)}

        self.initial_node = None
        self.final_node = None

        self.initialise_graph()

    def validate_codon_restrictions(self, codon_restrictions: Dict[int, CodonRestriction]) -> Dict[int, List[str]]:
        """
        Check the inputted codon restrictions make sense!
        """
        codon_restrictions = codon_restrictions or {}
        normalised = {}

        for pos, codon_restriction in codon_restrictions.items():
            if pos < 1 or pos > len(self.aa_seq):
                raise ValueError(f"Restricted position {pos} is out of range.")

            if isinstance(codon_restriction, str):
                codons = [codon_restriction]
            else:
                codons = list(codon_restriction)

            if len(codons) == 0:
                raise ValueError(f"Codon restriction at position {pos} cannot be empty.")

            codons = [codon.upper() for codon in codons]

            aa = self.aa_seq[pos - 1]

            if pos in self.codon_restrictions:
                allowed_codons = [self.ct.normalise_codon(codon) for codon in self.codon_restrictions[pos]]
            else:
                allowed_codons = self.ct.aa_to_codons[aa]

            for codon in codons:
                if codon not in allowed_codons:
                    raise ValueError(f'Codon {codon} is not valid for amino acid {aa} at position {pos}.')

            normalised[pos] = codons

        return normalised

    def initialise_graph(self) -> None:
        """
        Initialise the codon graph.
        """
        left_context_node = ContextNode(self.context_l)
        right_context_node = ContextNode(self.context_r)
        end_node = EndNode()

        for ix, aa in enumerate(self.aa_seq):
            pos = ix + 1

            if pos in self.codon_restrictions:
                codons = self.codon_restrictions[pos]
            else:
                codons = self.ct.aa_to_codons[aa]

            node = CodonNode(pos, aa, codons)
            self.add_codon_node(node)

        # Left context -> first codon node
        left_context_node.transitions = {
            left_context_node.sequence: self.codon_nodes[0]
        }
        self.codon_nodes[0].parents.add(
            (left_context_node, left_context_node.sequence)
        )

        # Codon node -> next codon node
        for i in range(1, len(self.codon_nodes)):
            previous = self.codon_nodes[i - 1]
            current = self.codon_nodes[i]

            for codon in previous.codons:
                previous.transitions[codon] = current
                current.parents.add((previous, codon))

        # Last codon node -> right context
        last_codon_node = self.codon_nodes[-1]
        for codon in last_codon_node.codons:
            last_codon_node.transitions[codon] = right_context_node
            right_context_node.parents.add((last_codon_node, codon))

        # Right context -> end
        right_context_node.transitions = {
            right_context_node.sequence: end_node
        }
        end_node.parents.add(
            (right_context_node, right_context_node.sequence)
        )

        self.left_context_node = left_context_node
        self.right_context_node = right_context_node
        self.end_node = end_node

        self.initial_node = left_context_node
        self.final_node = end_node

    @property
    def nodes(self) -> List[Node]:
        """
        All nodes in the graph, including context and end nodes.
        """
        return [
            self.left_context_node,
            *self.codon_nodes,
            self.right_context_node,
            self.end_node,
        ]

    def add_codon_node(self, node: CodonNode) -> None:
        """
        Add a codon node and update codon-node indexes.
        """
        self.codon_nodes.append(node)
        self.codon_nodes_by_pos.setdefault(node.pos, []).append(node)

    def remove_codon_node(self, node: CodonNode) -> None:
        """
        Remove a codon node and update codon-node indexes.
        """
        self.codon_nodes.remove(node)
        self.codon_nodes_by_pos[node.pos].remove(node)

    def view(self) -> "CodonGraphView":
        """
        Return a constrained view over this graph.
        """
        return CodonGraphView(self)


class CodonGraphView:
    """
    View of a codon graph. The view allows optional temporary constraints to be added without
    affecting the underlying codon graph. It is on this object that most operations take place.
    """

    def __init__(self,
                 graph: CodonGraph,
                 ) -> None:
        """
        Constructor for the CodonGraphView

        Parameters
        ----------
        graph
            The underlying codon graph.
        """
        self.graph = graph
        self.pinned_codons: Dict[int, List[str]] = {}

        self.valid_paths_by_choice = {}
        self.n_valid_sequences = None
        self.samplers = {}

        self.compile()

    @property
    def aa_seq(self):
        """
        The amino acid sequence on the underlying graph.

        Returns
        -------
        The aa seq.
        """
        return self.graph.aa_seq

    def __getitem__(self, index: int) -> str:
        """
        Return the valid sequence at a given index.

        Parameters
        ----------
        index
            Zero-based sequence index.

        Returns
        -------
        str
            The indexed valid DNA sequence.
        """

        if index < 0 or index >= self.n_valid_sequences:
            raise IndexError(f"Sequence index {index} out of range for {self.n_valid_sequences} valid sequences.")

        node = self.graph.initial_node
        sequence = []

        while node is not self.graph.final_node:
            choice_counts = self.valid_paths_by_choice[node]

            if isinstance(node, CodonNode):
                remaining = index

                for choice, count in choice_counts.items():
                    if remaining < count:
                        emitted = choice
                        break
                    remaining -= count
                else:
                    raise RuntimeError("Failed to resolve sequence index.")

                sequence.append(emitted)
                index = remaining

            else:
                # Context nodes only have one valid outgoing choice.
                emitted = next(iter(choice_counts))

            node = node.transitions[emitted]

        return ''.join(sequence)

    def __len__(self):
        """
        The number of valid sequences in this graph.

        Returns
        -------
        The number of valid sequences in this graph.
        """
        return self.n_valid_sequences

    def __iter__(self) -> Generator[str, None, None]:
        """
        Iterate over all valid sequences in this graph view.

        Yields
        ----------
        All valid sequences in the graph, in order.
        """
        yield from self.enumerate()

    def pin_codons(self, pinned_codons: Dict[int, CodonRestriction]) -> None:
        """
        Pin (temporarily fix) a codon in this codon graph view

        Parameters
        ----------
        pinned_codons
            A dict specifying which codons to pin, by pos: codon.
        """
        pinned_codons = self.graph.validate_codon_restrictions(pinned_codons)
        self.pinned_codons.update(pinned_codons)
        self.compile()

    def unpin_codons(self, positions: Sequence[int]) -> None:
        """
        Unpin codon nodes by pos.

        Parameters
        ----------
        positions
            A list of positions to unpin.
        """
        for pos in positions:
            if pos not in self.graph.codon_nodes_by_pos:
                raise ValueError(f'Pinned codon position {pos} is out of range.')

            self.pinned_codons.pop(pos, None)

        self.compile()

    def clear_pins(self) -> None:
        """
        Remove all codon pins from this graph view
        """
        self.pinned_codons.clear()
        self.compile()

    def contains(self, seq: str) -> bool:
        """
        Check whether a DNA sequence is contained in this view.

        Parameters
        ----------
        seq
            The sequence to check

        Returns
        -------
        True if and only if the sequence is contained in this sequence space.
        """
        seq = seq.upper()

        if len(seq) != len(self.graph.aa_seq) * 3:
            return False

        current_node = self.graph.initial_node.transitions[self.graph.initial_node.sequence]

        while current_node is not self.graph.right_context_node:
            pos = current_node.pos
            codon = seq[(pos - 1) * 3: pos * 3]

            if pos in self.pinned_codons:
                if codon not in self.pinned_codons[pos]:
                    return False

            elif codon not in current_node.codons:
                return False

            current_node = current_node.transitions[codon]

        return True

    def sample(self) -> str:
        """
        Sample a DNA sequence from this graph view.
        """
        node = self.graph.initial_node
        sequence = []

        while node is not self.graph.final_node:
            if not node.transitions:
                raise RuntimeError(f"Reached non-final node {node.id} with no outgoing transitions.")

            if isinstance(node, CodonNode):
                emitted = self.samplers[node].sample()
                sequence.append(emitted)

            else:
                emitted = node.sequence

            node = node.transitions[emitted]

        return ''.join(sequence)

    def enumerate(self, include_context: bool = False) -> Generator[str, None, None]:
        """
        Enumerate all valid sequences in this view.

        Parameters
        ----------
        include_context
            Whether to include left and right context sequences.

        Yields
        ------
        str
            A valid DNA sequence.
        """
        for index in range(self.n_valid_sequences):
            yield self.__getitem__(index)

    def copy(self) -> "CodonGraphView":
        """
        Copy this view and all its constraints.

        Returns
        -------
        A copy of the view.
        """
        view = self.graph.view()
        view.pin_codons(self.pinned_codons.copy())
        return view

    def compile(self) -> None:
        """
        Calculate all graph properties that are derived from its structure plus additional
        temporary constraints (pins). Remember to do this after editing the graph!
        """
        # Calculate descendant counts!
        self._update_descendant_counts()

        # Update the samplers!
        self._update_samplers()

    def _update_samplers(self) -> None:
        """
        Make samplers for each codon node given the view's constraints. The probabilities
        are calculated by combining descendant counts and the base probabilities.
        """
        samplers = {}

        for node, choice_counts in self.valid_paths_by_choice.items():
            if not isinstance(node, CodonNode):
                continue

            probability_by_codon = dict(zip(node.codons, node.probabilities))

            codons = []
            weights = []

            for codon, descendant_count in choice_counts.items():
                codons.append(codon)
                weights.append(probability_by_codon[codon] * descendant_count)

            if codons:
                samplers[node] = Sampler(codons, weights)

        self.samplers = samplers

    def _choices_for_node(self, node: CodonNode) -> List[str]:
        """
        Return codon choices available to this node in this view.
        """
        if node.pos in self.pinned_codons:
            return self.pinned_codons[node.pos]

        return node.codons

    def _update_descendant_counts(self) -> None:
        """
        Calculate valid path counts for each outgoing transition.
        """
        valid_paths_by_choice = {}

        right_choice = self.graph.right_context_node.sequence
        valid_paths_by_choice[self.graph.right_context_node] = {right_choice: 1}
        next_counts = {self.graph.right_context_node: 1}

        for pos in range(len(self.graph.aa_seq), 0, -1):
            current_counts = {}

            for node in self.graph.codon_nodes_by_pos[pos]:
                choice_counts = {}
                total = 0

                for choice in self._choices_for_node(node):
                    child = node.transitions[choice]
                    count = next_counts.get(child, 0)

                    if count:
                        choice_counts[choice] = count
                        total += count

                valid_paths_by_choice[node] = choice_counts
                current_counts[node] = total

            next_counts = current_counts

        left_choice = self.graph.left_context_node.sequence
        left_child = self.graph.left_context_node.transitions[left_choice]
        total = next_counts.get(left_child, 0)

        valid_paths_by_choice[self.graph.left_context_node] = {left_choice: total}

        self.valid_paths_by_choice = valid_paths_by_choice
        self.n_valid_sequences = total
