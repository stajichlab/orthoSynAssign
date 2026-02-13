from __future__ import annotations

from collections import defaultdict, deque
from itertools import combinations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ._gene import Gene, Genome, Protein


class Orthogroup:
    """Represents a group of orthologous genes."""

    __slots__ = ("og_id", "_genes")

    def __init__(self, og_id: str):
        self.og_id = og_id
        self._genes: list[Gene] = []

    def __len__(self) -> int:
        return len(self._genes)

    def __getitem__(self, index: int) -> Gene:
        return self._genes[index]

    def __iter__(self):
        return iter(self._genes)

    def add_gene(self, gene_obj: Gene):
        if gene_obj not in self._genes:
            gene_obj.orthogroup = self
            self._genes.append(gene_obj)

    def add_protein(self, protein_obj: Protein):
        self.add_gene(protein_obj.gene)

    def get_samples(self):
        return {gene.genome for gene in self}

    def has_paralogs(self):
        """Checks if any genome is represented more than once."""
        from collections import Counter

        counts = Counter(gene.genome for gene in self)
        return any(count > 1 for count in counts.values())

    def get_refined_sogs(self, window_size, ratio_threshold, skip_single_ortholog=True):
        """
        Refines the Orthogroup into SOGs based on synteny.
        """
        # Quick exit for non-paralogous groups
        if skip_single_ortholog and not self.has_paralogs():
            return [self._as_sog_dict(self._genes)]

        # Find the syntenic 'edges' using pairwise comparisons
        refined_pairs = self._perform_pairwise_comparisons(window_size, ratio_threshold)

        # If no pairs pass, the OG is essentially unsupported/fragmented
        if not refined_pairs:
            return []

        # Cluster the pairs into SOGs (via BFS/Connected Components)
        return _consolidate_into_sogs(refined_pairs)

    def _perform_pairwise_comparisons(self, window_size, ratio_threshold):
        """Internal helper to organize genes by genome and run comparisons."""
        genes_by_genome = defaultdict(list)
        for gene in self:
            genes_by_genome[gene.genome].append(gene)

        all_refined_pairs = []
        # Combinations ensures we only compare Genome A to Genome B once
        for genome_a, genome_b in combinations(genes_by_genome.keys(), 2):
            pairs = compare_gene_sets(genes_by_genome[genome_a], genes_by_genome[genome_b], window_size, ratio_threshold)
            all_refined_pairs.extend(pairs)
        return all_refined_pairs

    def _as_sog_dict(self, gene_list: list[Gene]) -> dict[str, list]:
        """Helper to format a list of genes into the final SOG output structure."""
        sog = defaultdict(list)
        for g in gene_list:
            sog[g.genome].append(g)
        return dict(sog)


def compare_gene_sets(
    genes_a: list[Gene], genes_b: list[Gene], window_size: int, ratio_threshold: float
) -> list[tuple[Gene, Gene]]:
    """
    Compare two sets of genes from different genomes within the same Orthogroup
    and return a list of syntenically supported orthologous pairs.
    """
    # Get the shared orthogroups once for this pair
    genome_a = genes_a[0].genome
    genome_b = genes_b[0].genome
    shared_ogs = get_shared_ogs(genome_a, genome_b)

    # Determine the primary and secondary genomes
    primary, secondary = (genes_a, genes_b) if len(genes_a) <= len(genes_b) else (genes_b, genes_a)

    # Pre-calculate all windows for the secondary set
    secondary_windows = {gene: gene.genome.get_window(gene, shared_ogs, window_size) for gene in secondary}

    refined_pairs = []

    # Compare windows
    for focal_gene in primary:
        win_p = focal_gene.genome.get_window(focal_gene, shared_ogs, window_size)

        # Find the best match
        best_candidate = None
        max_r = -1.0

        for candidate, win_s in secondary_windows.items():
            ratio = _calculate_synteny_ratio(win_p, win_s)

            if ratio >= ratio_threshold and ratio > max_r:
                max_r = ratio
                best_candidate = candidate

        if best_candidate:
            refined_pairs.append((focal_gene, best_candidate))

    return refined_pairs


def get_shared_ogs(genome_a: Genome, genome_b: Genome) -> set[str]:
    """
    Returns a set of OG IDs that are present in both genomes.
    This acts as the 'common language' for the comparison.
    """
    ogs_a = {g.orthogroup.og_id for g in genome_a if g.orthogroup}
    ogs_b = {g.orthogroup.og_id for g in genome_b if g.orthogroup}

    return ogs_a.intersection(ogs_b)


def _calculate_synteny_ratio(win_a: list[str], win_b: list[str]) -> float:
    """
    Calculates the 1-to-1 synteny match ratio between two dynamic windows.
    """
    if not win_a or not win_b:
        return 0.0

    # Count 1-to-1 matches
    matches = 0
    remaining_b = list(win_b)

    for og_id in win_a:
        if og_id in remaining_b:
            matches += 1
            remaining_b.remove(og_id)

    # Calculate the ratio of matches to total genes
    denominator = max(len(win_a), len(win_b))

    return matches / denominator


def _consolidate_into_sogs(pairs: list[tuple[Gene, Gene]]):
    """
    Standard Connected Components algorithm to cluster genes.
    Expects a list of gene pairs: [(gene1, gene2), (gene3, gene4), ...]
    """
    # Build an adjacency list (the graph)
    adj = defaultdict(list)
    nodes = set()
    for u, v in pairs:
        adj[u].append(v)
        adj[v].append(u)
        nodes.update([u, v])

    visited = set()
    sog_clusters = []

    # Traverse the graph to find connected components
    for node in nodes:
        if node not in visited:
            # Start a new SOG
            component = []
            queue = deque([node])
            visited.add(node)

            while queue:
                curr = queue.popleft()
                component.append(curr)

                for neighbor in adj[curr]:
                    if neighbor not in visited:
                        visited.add(neighbor)
                        queue.append(neighbor)

            # Format the cluster
            sog_clusters.append(_format_as_sog_dict(component))

    return sog_clusters


def _format_as_sog_dict(gene_list: list[Gene]) -> dict[Genome, str]:
    """Formats a list of genes into {Genome: [Genes]}."""
    sog = defaultdict(list)
    for g in gene_list:
        sog[g.genome].append(g)
    return dict(sog)


def align_sog_dict(sog_dict: dict[Gene, list[Gene]]) -> dict[Gene, list[Gene]]:
    """
    Aligns the neighborhood lists within the dictionary so that
    the focal gene (the key) is at the same index in every value.
    """
    gap_char = None
    if not sog_dict:
        return {}

    # Calculate the 'pivot' (the max index of the key in its value list)
    # Store the current index of each key for efficiency
    offsets = {}
    max_prefix_len = 0

    for focal_gene, neighborhood in sog_dict.items():
        try:
            # Where is the focal gene in its own neighborhood?
            idx = neighborhood.index(focal_gene)
        except ValueError:
            # Fallback if focal_gene isn't in the list
            idx = 0

        offsets[focal_gene] = idx
        if idx > max_prefix_len:
            max_prefix_len = idx

    # Build the aligned dictionary
    aligned_dict = {}

    # Calculate global max length to equalize the 'tails' later
    # First pass: calculate how long each list will be after front padding
    temp_lengths = []
    for focal_gene, neighborhood in sog_dict.items():
        front_pad_size = max_prefix_len - offsets[focal_gene]
        temp_lengths.append(front_pad_size + len(neighborhood))

    max_total_len = max(temp_lengths)

    # Apply padding
    for focal_gene, neighborhood in sog_dict.items():
        # Front padding: shifts the focal gene to the 'pivot' column
        front_pad = [gap_char] * (max_prefix_len - offsets[focal_gene])

        # Back padding: ensures all lists have the same total length
        current_padded_list = front_pad + list(neighborhood)
        back_pad = [gap_char] * (max_total_len - len(current_padded_list))

        aligned_dict[focal_gene] = current_padded_list + back_pad

    return aligned_dict
