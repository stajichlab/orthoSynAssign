from __future__ import annotations

from collections import defaultdict
from itertools import combinations

import numpy as np

from ._dsu import DSU
from .gene import Genome
from .orthogroup import Orthogroup


class SyntenyEngine:
    __slots__ = ("og_arrays", "seq_arrays", "og_members", "og_rep_maps", "shared_matrix")

    def __init__(self, genomes: list[Genome], orthogroups: list[Orthogroup]):
        # String -> Int for internal math; Int -> String for final reporting
        og_str_to_int = {og.id: i for i, og in enumerate(orthogroups)}

        genome_to_idx = {genome.name: i for i, genome in enumerate(genomes)}

        # Initialize data containers
        self.og_arrays: list[np.ndarray] = []  # List of int-represented gene arrays (og id, one array per genome)
        self.seq_arrays: list[np.ndarray] = []  # List of int-represented gene arrays (scaffold boundaries, one array per genome)

        # Populate data containers
        for genome in genomes:
            # Map seqid strings to local integers (e.g., "Chr1" -> 0)
            seqid_map = {s: i for i, s in enumerate(dict.fromkeys(gene.seqid for gene in genome._genes))}

            # [1, 2, 3,...] for OGs in a certain genome
            self.og_arrays.append(
                np.array([og_str_to_int.get(g.og.id, -1) if g.og else -1 for g in genome._genes], dtype=np.int32)
            )
            # [1, 1, 1,...] for scaffolds in a certain genome
            self.seq_arrays.append(np.array([seqid_map[g.seqid] for g in genome._genes], dtype=np.int32))

        # Parallel Mapping Lists
        self.og_members: list[list[tuple[int, str]]] = [[] for _ in range(len(orthogroups))]
        # member_to_backbone[og_idx][gene_id] = (genome_idx, backbone_idx)
        self.og_rep_maps: list[dict[tuple[int, str], tuple[int, int]]] = [{} for _ in range(len(orthogroups))]

        for og_idx, og in enumerate(orthogroups):
            for gene_obj in og:
                genome_idx = genome_to_idx[gene_obj.genome.name]
                # The physical anchor is always the representative's index in _genes
                rep_gene_idx = gene_obj.representative.index

                gene_id = gene_obj.id
                self.og_members[og_idx].append((genome_idx, gene_id))
                self.og_rep_maps[og_idx][(genome_idx, gene_id)] = (genome_idx, rep_gene_idx)

        # Pre-calculate the shared OG matrix
        self.shared_matrix = self._build_shared_matrix()

    def refine(self, og_idx: int, window_size: int, ratio_threshold: float) -> list[list[tuple[int, str]]]:
        """
        Coordinates the refinement of a single Orthogroup.
        Used by multiprocessing workers.
        """
        # Get the members for this OG
        member_ids = self.og_members[og_idx]
        if not member_ids:
            return []

        # Get the representative gene for this OG
        og_rep_map = self.og_rep_maps[og_idx]

        # Group isoforms by their shared backbone anchor
        # anchor_to_ids[(g_idx, b_idx)] = ["iso1", "iso2"]
        anchor_to_ids: defaultdict[tuple[int, int], list[tuple[int, str]]] = defaultdict(list)
        for m_id in member_ids:
            anchor_to_ids[og_rep_map[m_id]].append(m_id)

        anchors = list(anchor_to_ids.keys())

        # Perform pairwise comparisons using internal integer logic
        refined_pairs = self._perform_pairwise_comparisons(anchors, window_size, ratio_threshold)
        clusters = cluster_refined_ogs(refined_pairs, anchors)

        # EXPAND: Turn clusters of anchors into isoform IDs
        results: list[list[tuple[int, str]]] = []
        for nodes in clusters:
            if len(nodes) <= 1:
                continue

            iso_ids = []
            for anchor in nodes:
                iso_ids.extend(anchor_to_ids[anchor])

            results.append(iso_ids)

        return results

    def _build_shared_matrix(self) -> list[list[set[int]]]:
        """Pre-calculates which OGs are shared between every pair of genomes."""
        num_genomes = len(self.og_arrays)
        # Convert each genome's array to a set once
        genome_sets = [set(arr[arr != -1]) for arr in self.og_arrays]

        # Build a 2D matrix of intersections
        return [[genome_sets[i] & genome_sets[j] for j in range(num_genomes)] for i in range(num_genomes)]

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
        shared_ogs = self.shared_matrix[p_idx][s_idx]

        # Pre-calculate windows for the second set
        secondary_windows = []
        s_og_array = self.og_arrays[s_idx]
        s_seq_array = self.seq_arrays[s_idx]
        for gene_idx in secondary:
            win_s = get_window(s_og_array, s_seq_array, gene_idx, shared_ogs, window_size)
            # This is your "if g.og.id in shared_ogs" logic, but with integers
            secondary_windows.append((gene_idx, win_s))

        refined_pairs = []

        # Compare
        p_og_array = self.og_arrays[p_idx]
        p_seq_array = self.seq_arrays[p_idx]
        for gene_idx in primary:
            win_p = get_window(p_og_array, p_seq_array, gene_idx, shared_ogs, window_size)

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


def get_window(og_array: np.ndarray, seq_array: np.ndarray, gene_idx: int, target_ogs: set[int], window_size: int) -> np.ndarray:
    """
    Retrieves a window of Orthogroup IDs, restricted by scaffold boundaries.
    """
    half_win = window_size // 2
    focal_seqid = seq_array[gene_idx]

    is_anchor = np.isin(og_array, list(target_ogs))
    is_same_scaffold = seq_array == focal_seqid

    valid_indices = np.where(is_anchor & is_same_scaffold)[0]

    pos = np.searchsorted(valid_indices, gene_idx)

    # Fetch raw slices
    left_side = valid_indices[max(0, pos - half_win) : pos]
    start_r = pos + 1 if (pos < len(valid_indices) and valid_indices[pos] == gene_idx) else pos
    right_side = valid_indices[start_r : start_r + half_win]

    neighbor_indices = np.concatenate([left_side, right_side])

    return og_array[neighbor_indices]


def calculate_synteny_ratio(win_a: np.ndarray, win_b: np.ndarray) -> float:
    """Calculates the 1-to-1 synteny match ratio between two dynamic windows.

    Specifically, this function computes the ratio of overlapping orthogroups present in both window sets. The overlap is
    determined by finding the minimum count for each shared Orthogroup ID across the two windows.

    Args:
        win_a (np.ndarray): An array of Orthogroup indices representing the first dynamic window.
        win_b (np.ndarray): An array of Orthogroup indices representing the second dynamic window.

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
