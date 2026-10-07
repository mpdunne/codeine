from typing import Dict, Iterable, Mapping, Sequence, Set, Tuple

from codeine.graph.factors import ChoiceFactor


def _dependency_graph(scopes: Iterable[Iterable[int]]) -> Dict[int, Set[int]]:
    """
    Map each position to the other positions with which it shares a factor.
    """
    neighbours = {}

    for scope in scopes:
        positions = tuple(scope)

        for pos in positions:
            neighbours.setdefault(pos, set())

        for ix, pos in enumerate(positions):
            neighbours[pos].update(positions[:ix])
            neighbours[pos].update(positions[ix + 1:])

    return neighbours


def _min_fill_order(scopes: Iterable[Iterable[int]]) -> Tuple[int, ...]:
    """
    Choose a diagram order that adds as few dependencies as possible.

    Removing a position connects its remaining neighbours. Prefer positions
    requiring the fewest new connections, then the fewest neighbours, then the
    lowest graph position. This keeps intermediate diagrams small in many cases.
    """
    neighbours = _dependency_graph(scopes)
    order = []

    while neighbours:
        def score(pos):
            adjacent = neighbours[pos]
            adjacent_tuple = tuple(adjacent)
            missing_edges = sum(
                neighbour_b not in neighbours[neighbour_a]
                for ix, neighbour_a in enumerate(adjacent_tuple)
                for neighbour_b in adjacent_tuple[ix + 1:]
            )
            return missing_edges, len(adjacent), pos

        pos = min(neighbours, key=score)
        adjacent = set(neighbours[pos])
        order.append(pos)

        for neighbour_a in adjacent:
            for neighbour_b in adjacent:
                if neighbour_a != neighbour_b:
                    neighbours[neighbour_a].add(neighbour_b)

        for neighbour in adjacent:
            neighbours[neighbour].discard(pos)

        del neighbours[pos]

    return tuple(order)


class _DecisionDiagrams:
    """
    Shared decision diagrams for individual factors.

    A node branches on one graph position, with one child per possible choice.
    FALSE and TRUE are rejecting and accepting leaves. Identical nodes are
    shared, and nodes whose choices all lead to the same result are removed.
    Separate factors retain separate roots: joining them all into one diagram
    would recreate the large product of constraint states we want to avoid.
    """

    FALSE = 0
    TRUE = 1

    def __init__(self, domains, variable_order=None):
        self.domains = {pos: tuple(choices) for pos, choices in domains.items()}
        if variable_order is None:
            variable_order = tuple(sorted(self.domains))
        ordered = tuple(variable_order) + tuple(pos for pos in sorted(self.domains) if pos not in variable_order)
        self.rank = {pos: ix for ix, pos in enumerate(ordered)}
        self.choice_ix = {
            pos: {choice: ix for ix, choice in enumerate(choices)}
            for pos, choices in self.domains.items()
        }
        self.nodes = [None, None]
        self.unique = {}
        self._apply_cache = {}
        self._not_cache = {self.FALSE: self.TRUE, self.TRUE: self.FALSE}
        self._restrict_cache = {}

    def var(self, node):
        return None if node < 2 else self.nodes[node][0]

    def children(self, node):
        return self.nodes[node][1]

    def make_node(self, pos, children):
        children = tuple(children)

        if all(child == children[0] for child in children[1:]):
            return children[0]

        key = (pos, children)
        node = self.unique.get(key)

        if node is None:
            node = len(self.nodes)
            self.unique[key] = node
            self.nodes.append(key)

        return node

    def negate(self, node):
        cached = self._not_cache.get(node)
        if cached is not None:
            return cached

        pos, children = self.nodes[node]
        result = self.make_node(pos, (self.negate(child) for child in children))
        self._not_cache[node] = result
        return result

    def conjunction(self, left, right):
        if left == self.FALSE or right == self.FALSE:
            return self.FALSE
        if left == self.TRUE:
            return right
        if right == self.TRUE or left == right:
            return left

        key = (left, right) if left < right else (right, left)
        cached = self._apply_cache.get(key)
        if cached is not None:
            return cached

        left_pos = self.var(left)
        right_pos = self.var(right)
        pos = left_pos if self.rank[left_pos] < self.rank[right_pos] else right_pos
        left_children = self.children(left) if left_pos == pos else None
        right_children = self.children(right) if right_pos == pos else None
        n_choices = len(self.domains[pos])

        result = self.make_node(
            pos,
            (
                self.conjunction(
                    left_children[ix] if left_children is not None else left,
                    right_children[ix] if right_children is not None else right,
                )
                for ix in range(n_choices)
            ),
        )
        self._apply_cache[key] = result
        return result

    def relation(self, reference_pos, compare_pos, allowed_by_reference_choice):
        """
        Build a decision diagram for one relation between two graph positions.
        """
        positions = tuple(sorted({reference_pos, compare_pos}, key=self.rank.__getitem__))

        def build(depth, assignments):
            if depth == len(positions):
                reference_choice = assignments[reference_pos]
                compare_choice = assignments[compare_pos]
                return self.TRUE if compare_choice in allowed_by_reference_choice.get(reference_choice, ()) \
                    else self.FALSE

            pos = positions[depth]
            children = []
            for choice in self.domains[pos]:
                assignments[pos] = choice
                children.append(build(depth + 1, assignments))
            assignments.pop(pos, None)
            return self.make_node(pos, children)

        return build(0, {})

    def allowed_factor(self, relations):
        """
        Build a factor that rejects only when every relation matches.
        """
        forbidden = self.TRUE

        for relation in relations:
            reference_pos, compare_pos, allowed = relation
            node = self.relation(reference_pos, compare_pos, allowed)

            forbidden = self.conjunction(forbidden, node)
            if forbidden == self.FALSE:
                return self.TRUE

        return self.negate(forbidden)

    def advance(self, node, pos, choice):
        """
        Condition a diagram on one graph-position choice.
        """
        if node < 2:
            return node

        key = (node, pos, choice)
        cached = self._restrict_cache.get(key)
        if cached is not None:
            return cached

        node_pos, children = self.nodes[node]
        if node_pos == pos:
            result = children[self.choice_ix[pos][choice]]
        elif self.rank[node_pos] > self.rank[pos]:
            result = node
        else:
            result = self.make_node(
                node_pos,
                (self.advance(child, pos, choice) for child in children),
            )

        self._restrict_cache[key] = result
        return result


class ComponentModelCounter:
    """
    Count valid graph assignments by solving independent groups of rules.

    Independent groups are AND branches: all must be satisfied, so their counts
    multiply. Choices within a connected group are OR branches: each choice
    contributes its number of valid completions, so their counts add.

    After each choice, factors are reduced and checked for new independent
    groups. Equivalent remaining problems share a cached count. Positions no
    longer used by any factor contribute their number of available choices.
    """

    def __init__(self, domains: Mapping[int, Sequence[str]], factors: Sequence[ChoiceFactor]) -> None:
        """
        Parameters
        ----------
        domains
            Possible choices at each graph position, including context nodes.
            Restricted choices can be supplied here to count a pinned space.
        factors
            Forbidden combinations produced by linked constraints.
        """
        factors = tuple(factors)
        self.manager = _DecisionDiagrams(domains, _min_fill_order(factor.scope for factor in factors))
        self.domains = self.manager.domains
        self._support_cache = {}
        self.n_calls = 0
        self.n_splits = 0

        # Each distinct set of remaining rules gets an ID. Cache its count and
        # the positions it still depends on, represented as a position bitmask.
        self._component_ids = {}
        self._component_roots = []
        self._component_supports = []
        self._component_counts = []
        self._component_splits = {}
        self._component_choices = {}
        self._component_transitions = {}

        self._initial_count = None

        if any(not choices for choices in self.domains.values()):
            initial_roots = (self.manager.FALSE,)
        else:
            initial_roots = tuple(
                self.manager.allowed_factor(factor.relations)
                for factor in factors
            )

        self._initial_component_id = self._get_component_id(initial_roots)
        all_variables = sum(1 << pos for pos in self.domains)
        self._initial_free_variables = all_variables & ~self._component_supports[self._initial_component_id]

    def _support(self, root):
        cached = self._support_cache.get(root)
        if cached is not None:
            return cached

        if root < 2:
            result = 0
        else:
            pos, children = self.manager.nodes[root]
            result = 1 << pos
            for child in children:
                result |= self._support(child)

        self._support_cache[root] = result
        return result

    def _canonical_roots(self, roots):
        return tuple(sorted(set(
            root
            for root in roots
            if root != self.manager.TRUE
        )))

    def _get_component_id(self, roots):
        roots = self._canonical_roots(roots)
        component_id = self._component_ids.get(roots)
        if component_id is not None:
            return component_id

        component_id = len(self._component_roots)
        support = 0
        for root in roots:
            support |= self._support(root)
        self._component_ids[roots] = component_id
        self._component_roots.append(roots)
        self._component_supports.append(support)
        self._component_counts.append(None)
        return component_id

    @staticmethod
    def _positions(mask):
        while mask:
            bit = mask & -mask
            yield bit.bit_length() - 1
            mask ^= bit

    def _free_count(self, variables):
        result = 1
        for pos in self._positions(variables):
            result *= len(self.domains[pos])
        return result

    def _components(self, roots, scopes):
        unseen = set(range(len(roots)))
        components = []

        while unseen:
            root_ix = unseen.pop()
            component = {root_ix}
            variables = scopes[root_ix]

            while True:
                adjacent = {
                    other_ix
                    for other_ix in unseen
                    if variables & scopes[other_ix]
                }
                if not adjacent:
                    break

                unseen.difference_update(adjacent)
                component.update(adjacent)
                for other_ix in adjacent:
                    variables |= scopes[other_ix]

            components.append(tuple(roots[ix] for ix in component))

        return components

    def _split_components(self, component_id):
        cached = self._component_splits.get(component_id)
        if cached is not None:
            return cached

        roots = self._component_roots[component_id]
        scopes = tuple(self._support(root) for root in roots)
        components = self._components(roots, scopes)
        result = tuple(self._get_component_id(component) for component in components)
        self._component_splits[component_id] = result
        return result

    def _choose_variable(self, component_id):
        cached = self._component_choices.get(component_id)
        if cached is not None:
            return cached

        occurrences = {}
        for root in self._component_roots[component_id]:
            for pos in self._positions(self._support(root)):
                occurrences[pos] = occurrences.get(pos, 0) + 1

        pos = max(
            occurrences,
            key=lambda candidate: (
                occurrences[candidate],
                -len(self.domains[candidate]),
                -candidate,
            ),
        )
        self._component_choices[component_id] = pos
        return pos

    def _advance_component(self, component_id, pos, choice):
        key = component_id, pos, choice
        cached = self._component_transitions.get(key)
        if cached is not None:
            return cached

        roots = self._component_roots[component_id]
        next_component_id = self._get_component_id(
            self.manager.advance(root, pos, choice)
            for root in roots
        )
        current_support = self._component_supports[component_id]
        next_support = self._component_supports[next_component_id]
        # Some rules become harmless after this choice. Their other positions
        # are now free and must still contribute their available choices.
        free_variables = current_support & ~(1 << pos) & ~next_support
        result = next_component_id, free_variables
        self._component_transitions[key] = result
        return result

    def _count_component(self, component_id):
        cached = self._component_counts[component_id]
        if cached is not None:
            return cached

        self.n_calls += 1
        roots = self._component_roots[component_id]

        if self.manager.FALSE in roots:
            result = 0
        elif not roots:
            result = 1
        else:
            components = self._split_components(component_id)
            if len(components) > 1:
                self.n_splits += 1
                result = 1
                for child_component_id in components:
                    result *= self._count_component(child_component_id)
            else:
                pos = self._choose_variable(component_id)
                result = 0
                for choice in self.domains[pos]:
                    child_component_id, free_variables = self._advance_component(
                        component_id, pos, choice
                    )
                    result += self._free_count(free_variables) * self._count_component(
                        child_component_id
                    )

        self._component_counts[component_id] = result
        return result

    def count(self) -> int:
        """
        Return the exact number of assignments satisfying all factors.
        """
        if self._initial_count is None:
            self._initial_count = self._free_count(self._initial_free_variables) * self._count_component(
                self._initial_component_id
            )
        return self._initial_count

    def contains(self, assignment: Mapping[int, str]) -> bool:
        """
        Check a complete assignment against the domains and every factor.

        Assignments must include all graph positions, including contexts, and
        no extra positions. This checks validity independently of sampling
        weights, without building new diagrams or changing counting caches.
        """
        if assignment.keys() != self.domains.keys():
            return False

        for pos, choice in assignment.items():
            if choice not in self.manager.choice_ix[pos]:
                return False

        for root in self._component_roots[self._initial_component_id]:
            while root >= 2:
                pos, children = self.manager.nodes[root]
                root = children[self.manager.choice_ix[pos][assignment[pos]]]

            if root == self.manager.FALSE:
                return False

        return True
