from typing import List, Sequence, Tuple, Union

from codeine.constraints.subpaths import SubPathConstraint
from codeine.graph.nodes import CodonNode
from codeine.motifs.restriction import RestrictionSite

from codeine.constraints.subpaths import SubPath


Motif = Union[str, RestrictionSite]
Motifs = Union[Motif, Sequence[Motif]]


class ForbiddenMotifConstraint(SubPathConstraint):
    """
    Prevent specified nucleotide motifs from occurring in generated sequences.

    Motifs may be provided as nucleotide strings or as ``RestrictionSite``
    objects. Restriction sites are expanded into their concrete recognition
    sequences during construction.
    """

    def __init__(self, forbidden_motifs: Motifs) -> None:
        """
        Parameters
        ----------
        forbidden_motifs
            One or more nucleotide motifs or restriction sites that must not occur.
        """
        super().__init__()

        if isinstance(forbidden_motifs, (str, RestrictionSite)):
            forbidden_motifs = [forbidden_motifs]

        self.motifs = tuple(forbidden_motifs)

        sequences = []

        for motif in forbidden_motifs:
            if isinstance(motif, RestrictionSite):
                sequences.extend(motif.motifs)
                continue

            if not isinstance(motif, str):
                raise TypeError('Forbidden motifs must be strings or RestrictionSite objects.')

            motif = motif.upper()

            if motif == '':
                raise ValueError('Forbidden motifs cannot be empty.')

            if not set(motif) <= set('ACGTU'):
                raise ValueError('Forbidden motifs must contain only A, C, G, T, or U.')

            sequences.append(motif)

        self.forbidden_sequences = tuple(sorted(set(sequences)))

    def _find_paths(self) -> Tuple[SubPath, ...]:
        """
        Find every concrete graph path that can generate a forbidden motif.

        Each returned SubPath records the emitted sequence, the graph steps
        required to produce it, and the offset at which the forbidden motif begins
        within the first emitted choice.

        Returns
        -------
        Tuple[SubPath, ...]
            All graph subpaths capable of producing one of the forbidden motifs.
        """
        paths = []

        sequences = [self.graph.tt.normalise_sequence(sequence) for sequence in self.forbidden_sequences]

        for sequence in sorted(set(sequences)):
            paths.extend(self._find_matching_subpaths(sequence))

        return tuple(paths)

    def _find_matching_subpaths(self, sequence: str) -> List[SubPath]:
        """
        Find all graph subpaths that can emit a given sequence.

        Parameters
        ----------
        graph
            The codon graph to search.
        sequence
            The sequence to search for.

        Returns
        -------
        Matching subpaths as a list of SubPath objects.
        """
        sequence = sequence.upper()
        graph = self.graph

        if len(sequence) == 0:
            raise ValueError('Sequence cannot be empty.')

        matches: List[SubPath] = []
        candidate_matches = []

        def add_match(partial_path, offset):
            steps = tuple(
                (node.pos, choice)
                for node, choice in partial_path
            )

            matches.append(
                SubPath(
                    sequence=sequence,
                    steps=steps,
                    offset=offset,
                )
            )

        # First, check which nodes we can start at.
        for node in graph.nodes:
            if node is graph.end_node:
                continue

            for choice in node.transitions:
                for offset in range(len(choice)):
                    choice_subsequence = choice[offset:]

                    if choice_subsequence.startswith(sequence):
                        # Bingo!
                        add_match(((node, choice),), offset)

                    elif sequence.startswith(choice_subsequence):
                        # Maybe bingo! Maygo!
                        candidate_matches.append((
                            ((node, choice),),
                            offset,
                            len(choice_subsequence),
                        ))

        def reinspect_candidate_matches(candidate_matches):
            reinspect = []

            for partial_path, offset, seen_length in candidate_matches:
                previous_node, previous_choice = partial_path[-1]
                node = previous_node.transitions[previous_choice]

                remaining_sequence = sequence[seen_length:]

                if remaining_sequence == '':
                    # Fantastic!
                    add_match(partial_path, offset)
                    continue

                if node is graph.end_node:
                    continue

                if isinstance(node, CodonNode):
                    choice_length = 3
                else:
                    choice_length = len(next(iter(node.transitions)))

                # Sneaky shortcut if we've crossed into the right context:
                if isinstance(node, CodonNode):
                    pos = node.pos

                    remaining_sequence_length = len(sequence) - seen_length
                    remaining_coding_length = 3 * (len(graph.aa_seq) - pos + 1)

                    if remaining_sequence_length > remaining_coding_length:
                        sequence_end = sequence[seen_length + remaining_coding_length:]

                        if not graph.context_r.startswith(sequence_end):
                            continue

                for choice in node.transitions:

                    if len(remaining_sequence) >= choice_length:

                        if remaining_sequence.startswith(choice):
                            # Keep going...
                            reinspect.append((
                                partial_path + ((node, choice),),
                                offset,
                                seen_length + choice_length,
                            ))

                        else:
                            # Hard luck this time.
                            continue

                    else:

                        if choice.startswith(remaining_sequence):
                            # Wahoo!
                            add_match(
                                partial_path + ((node, choice),),
                                offset,
                            )

                        else:
                            # Hard luck this time.
                            continue

            return reinspect

        while candidate_matches:
            candidate_matches = reinspect_candidate_matches(candidate_matches)

        return matches
