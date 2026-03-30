from __future__ import annotations

from collections import defaultdict
from itertools import combinations
from typing import TYPE_CHECKING

import numpy as np

from ._dsu import DSU

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from .gene import Genome
    from .orthogroup import Orthogroup


class SyntenyEngine:
    __slots__ = ("og_arrays", "seq_arrays", "orthogroups", "shared_og_matrix")

    def __init__(self, genomes: list[Genome], orthogroups: list[Orthogroup]):
        # String -> Int for internal math; Int -> String for final reporting
        og_str_to_int = {og.id: i for i, og in enumerate(orthogroups)}

        # Initialize data containers
        # List of int-represented arrays for og id, one array per genome)
        self.og_arrays: list[NDArray[np.int32]] = []
        # List of int-represented arrays for scaffold boundaries, one array per genome)
        self.seq_arrays: list[NDArray[np.int16]] = []
        self.orthogroups: list[list[tuple[int, str]]] = [[] for _ in range(len(orthogroups))]

        # Populate data containers
        for genome_idx, genome in enumerate(genomes):
            # Map seqid strings to local integers (e.g., "Chr1" -> 0)
            seqid_map: dict[str, int] = {}
            next_seqid_int = 0

            og_arr = []
            seqid_arr = []

            for gene_idx, gene in enumerate(genome._genes):
                seqid = gene.seqid
                if seqid not in seqid_map:
                    seqid_map[seqid] = next_seqid_int
                    next_seqid_int += 1
                # [1, 1, 1,...] for scaffolds in a certain genome
                seqid_arr.append(seqid_map[seqid])

                og_int = og_str_to_int.get(gene.og.id, -1) if gene.og else -1
                # [1, 2, 3,...] for OGs in a certain genome
                og_arr.append(og_int)

                # Populate the orthogroups for the engine's 'refine' task
                if og_int != -1:
                    self.orthogroups[og_int].append((genome_idx, gene_idx))

            # Convert to array
            self.og_arrays.append(np.array(og_arr, dtype=np.int32))
            self.seq_arrays.append(np.array(seqid_arr, dtype=np.int16))

        # Pre-calculate the shared OG matrix
        self.shared_og_matrix = self._build_shared_og_matrix()

    def refine(self, og_idx: int, window_size: int, ratio_threshold: float) -> list[list[tuple[int, str]]]:
        """
        Coordinates the refinement of a single Orthogroup using physical anchors.

        Returns:
            list[list[tuple[int, int]]]: A list of clusters, where each cluster is a
            list of (genome_idx, gene_idx) physical anchors.
        """
        # Get the genes for this OG
        genes = self.orthogroups[og_idx]
        if not genes:
            return []

        # Perform pairwise comparisons using internal integer logic
        refined_pairs = self._perform_pairwise_comparisons(genes, window_size, ratio_threshold)
        clusters = cluster_refined_ogs(refined_pairs, genes)

        return [cluster for cluster in clusters if len(cluster) > 1]

    def _build_shared_og_matrix(self) -> list[list[NDArray[np.int32]]]:
        """Pre-calculates which OGs are shared between every pair of genomes."""
        num_genomes = len(self.og_arrays)
        # Convert each genome's array to a set once
        genome_sets = [set(arr[arr != -1]) for arr in self.og_arrays]

        matrix = [[None] * num_genomes for _ in range(num_genomes)]

        for i in range(num_genomes):
            for j in range(i, num_genomes):
                # Calculate intersection using set logic
                intersection_set = genome_sets[i] & genome_sets[j]

                # CONVERT TO ARRAY: Sort it to speed up np.isin later
                # np.unique returns a sorted array which np.isin loves
                intersect_arr = np.array(list(intersection_set), dtype=np.int32)
                intersect_arr.sort()

                matrix[i][j] = intersect_arr
                matrix[j][i] = intersect_arr

        return matrix

    def _perform_pairwise_comparisons(self, members: list[tuple[int, int]], window_size: int, ratio_threshold: float):
        """Internal logic using the shared_matrix and og_arrays."""
        refined_pairs = []

        # Group members by genome to replicate your 'primary vs secondary' logic

        genes_by_genome = defaultdict(list)
        for genome_idx, gene_idx in members:
            genes_by_genome[genome_idx].append(gene_idx)

        genome_indices = list(genes_by_genome.keys())

        # All-vs-All Genome Pairs within this Orthogroup
        for genome_a, genome_b in combinations(genome_indices, 2):
            # Use your optimized 'compare_gene_pairs' logic here
            pairs = self._compare_gene_pairs(genome_a, genome_b, genes_by_genome, window_size, ratio_threshold)
            refined_pairs.extend(pairs)

        return refined_pairs

    def _compare_gene_pairs(
        self,
        genome_a: int,
        genome_b: int,
        genes_by_genomes: dict[int, list[int]],
        window_size: int,
        ratio_threshold: float,
    ) -> list[tuple[tuple[int, int], tuple[int, int]]]:
        genes_a, genes_b = genes_by_genomes[genome_a], genes_by_genomes[genome_b]

        # Compare lengths and unpack into fixed roles
        if len(genes_a) <= len(genes_b):
            primary, p_idx = genes_a, genome_a
            secondary, s_idx = genes_b, genome_b
        else:
            primary, p_idx = genes_b, genome_b
            secondary, s_idx = genes_a, genome_a

        # Get the shared OG set for this genome pair (from our pre-calculated matrix)
        shared_ogs = self.shared_og_matrix[p_idx][s_idx]

        # Pre-calculate windows for the second set
        secondary_windows = []
        s_og_array = self.og_arrays[s_idx]
        s_seq_array = self.seq_arrays[s_idx]
        s_og_array_masked = np.isin(s_og_array, shared_ogs)
        for gene_idx in secondary:
            win_s = s_og_array[get_window(s_og_array_masked, s_seq_array, gene_idx, window_size)]
            # This is your "if g.og.id in shared_ogs" logic, but with integers
            secondary_windows.append((gene_idx, win_s))

        refined_pairs = []

        # Compare
        p_og_array = self.og_arrays[p_idx]
        p_seq_array = self.seq_arrays[p_idx]
        p_og_array_masked = np.isin(p_og_array, shared_ogs)
        for gene_idx in primary:
            win_p = p_og_array[get_window(p_og_array_masked, p_seq_array, gene_idx, window_size)]

            best_candidate = -1
            max_r = -1.0

            for candidate_idx, win_s in secondary_windows:
                ratio = calculate_synteny_ratio(win_p, win_s)

                if ratio >= ratio_threshold and ratio > max_r:
                    max_r = ratio
                    best_candidate = candidate_idx

            if best_candidate != -1:
                p_coord = (p_idx, gene_idx)
                s_coord = (s_idx, best_candidate)
                refined_pairs.append((p_coord, s_coord))

        return refined_pairs


def get_window(og_masked_array: NDArray[np.bool_], seq_array: NDArray[np.int16], gene_idx: int, window_size: int) -> NDArray:
    """
    Retrieves a window of indices, restricted by scaffold boundaries.
    """
    half_win = window_size // 2
    focal_seqid = seq_array[gene_idx]

    is_same_scaffold = seq_array == focal_seqid

    valid_indices = np.where(og_masked_array & is_same_scaffold)[0]

    pos = np.searchsorted(valid_indices, gene_idx)

    # Fetch raw slices
    left_side = valid_indices[max(0, pos - half_win) : pos]
    start_r = pos + 1 if (pos < len(valid_indices) and valid_indices[pos] == gene_idx) else pos
    right_side = valid_indices[start_r : start_r + half_win]

    return np.concatenate([left_side, right_side])


def calculate_synteny_ratio(win_a: NDArray[np.int32], win_b: NDArray[np.int32]) -> float:
    """Calculates the 1-to-1 synteny match ratio between two dynamic windows.

    Specifically, this function computes the ratio of overlapping orthogroups present in both window sets. The overlap is
    determined by finding the minimum count for each shared Orthogroup ID across the two windows.

    Args:
        win_a (NDArray): An array of Orthogroup indices representing the first dynamic window.
        win_b (NDArray): An array of Orthogroup indices representing the second dynamic window.

    Returns:
        float: The synteny ratio, calculated as the number of overlapping orthogroups divided by the length of the longer window.
        If either window is empty, the function returns 0.0.
    """
    # Quick exit for empty windows (e.g., end of scaffold)
    len_a = win_a.size
    len_b = win_b.size
    if len_a == 0 or len_b == 0:
        return 0.0

    # Get unique OGs and their frequencies
    # For speed in the engine, assume -1 was already filtered out before calling this
    u_a, counts_a = np.unique(win_a, return_counts=True)
    u_b, counts_b = np.unique(win_b, return_counts=True)

    # Convert to dictionaries for O(1) intersection lookup
    # dict(zip(...)) is very efficient for small window sizes (e.g., 10-20 genes)
    dict_a = dict(zip(u_a, counts_a))
    dict_b = dict(zip(u_b, counts_b))

    # Calculate matches (min count of shared IDs)
    # Iterate over the smaller unique set to minimize lookups
    matches = 0
    small, large = (dict_a, dict_b) if len(dict_a) < len(dict_b) else (dict_b, dict_a)
    for og_id, count in small.items():
        if og_id in large:
            matches += min(count, large[og_id])

    # Normalize by the longer window length
    return matches / max(len_a, len_b)


def cluster_refined_ogs(
    pairs: list[tuple[tuple[int, int], tuple[int, int]]], all_genes: list[tuple[int, int]] = None
) -> list[list[tuple[int, int]]]:
    """
    Consolidates syntenic coordinate pairs into clusters using DSU.

    Args:
        pairs: List of ((g1, i1), (g2, i2)) representing syntenic links.
        all_genes: The full list of (genome_idx, gene_idx) in the Orthogroup.
                   Essential for including singletons that have no pairs.
    """
    # Collect all genes that need to be in the graph
    # If all_genes is provided, use it; otherwise, only use genes found in pairs.
    if all_genes is not None:
        nodes = all_genes
    else:
        nodes = list({coord for pair in pairs for coord in pair})

    if not nodes:
        return []

    # Map Genes to 0...N-1 indices for the DSU
    gene_to_id = {coord: i for i, coord in enumerate(nodes)}

    # Initialize your external DSU class
    dsu = DSU(len(nodes))

    # Process syntenic edges
    for u, v in pairs:
        dsu.union(gene_to_id[u], gene_to_id[v])

    # Group genes by their root parent
    clusters: defaultdict[int, list[tuple[int, int]]] = defaultdict(list)
    for i, coord in enumerate(nodes):
        root = dsu.find(i)
        clusters[root].append(coord)

    # Sort the members WITHIN each cluster
    # This ensures (Genome 0, Index 10) always comes before (Genome 1, Index 5)
    for root in clusters:
        clusters[root].sort()

    # Sort the clusters THEMSELVES
    # We sort the list of lists based on the first element of each sub-list
    sorted_clusters = sorted(clusters.values(), key=lambda x: x[0])

    return sorted_clusters
