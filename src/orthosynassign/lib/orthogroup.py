from __future__ import annotations

from collections import Counter, defaultdict, deque
from itertools import combinations
from typing import TYPE_CHECKING, Any, Iterator

if TYPE_CHECKING:
    from .gene import Gene, Genome


class Orthogroup:
    """Represents a group of orthologous genes.

    Attributes:
        id: The unique identifier for the orthogroup.
        _genes: A private attribute to store the genes in the orthogroup.
    """

    __slots__ = ("id", "_genes")

    def __init__(self, og_id: str | None = None, genes: list[Gene] | None = None) -> None:
        """Initialize an Orthogroup.

        Args:
            og_id (str): The unique identifier for the orthogroup.
            genes (list[Gene]): The genes in the orthogroup.
        """
        self.id = og_id
        self._genes: list[Gene] = []
        if genes:
            for gene in genes:
                self.add_gene(gene)

    def __repr__(self) -> str:
        """Return a string representation of the Orthogroup.

        Returns:
            str: A string in the format "[{id} | with {len(self._genes)} genes]".
        """
        og_id = self.id if self.id else "Unnamed orthogroup"
        return f"[{og_id} | with {len(self._genes)} genes]"

    def __len__(self) -> int:
        """Return the number of genes in the Orthogroup.

        Returns:
            int: The number of genes.
        """
        return len(self._genes)

    def __getitem__(self, index: int) -> Gene:
        """Retrieve a gene by its index in the Orthogroup.

        Args:
            index (int): The index of the gene to retrieve.

        Returns:
            Gene: The gene at the specified index.
        """
        return self._genes[index]

    def __contains__(self, item: Gene) -> bool:
        """Allows for 'gene in orthogroup' syntax."""
        return item in self._genes

    def __iter__(self) -> Iterator[Gene]:
        """Return an iterator over the genes in the Orthogroup.

        Returns:
            Iterator[Gene]: An iterator over the genes.
        """
        return iter(self._genes)

    def __getstate__(self) -> dict[str, Any]:
        """Get the state of the object for pickling.

        Returns:
            dict[str, Any]: The state of the object.
        """
        return {slot: getattr(self, slot) for slot in self.__slots__ if hasattr(self, slot)}

    def __setstate__(self, state: dict[str, Any]) -> None:
        """Restore the state of the object from pickling.

        Args:
            state (dict[str, Any]): The state to restore.
        """
        for slot, value in state.items():
            setattr(self, slot, value)

    def add_gene(self, gene: Gene) -> None:
        """Add a gene to the Orthogroup.

        If the gene is not already present in the Orthogroup, it will be added and its og attribute set to this Orthogroup.

        Args:
            gene (Gene): The gene to add.
        """
        if gene not in self._genes:
            gene.og = self
            self._genes.append(gene)

    def refine(self, window_size: int, ratio_threshold: float) -> list[SOG]:
        """Refine the orthologous gene groups (SOGs) within the orthogroup by finding syntenic 'edges' using pairwise comparisons.

        Args:
            window_size (int): The size of the window used for comparison.
            ratio_threshold (float): The threshold ratio used to determine if a pair of genes is considered syntenic.

        Returns:
            list[SOG]: A list of SOG (refined orthogroup) objects.
        """
        # Find the syntenic 'edges' using pairwise comparisons
        refined_pairs = _perform_pairwise_comparisons(self, window_size, ratio_threshold)

        # If no pairs pass, the OG is essentially unsupported/fragmented
        if not refined_pairs:
            return []

        # Cluster the pairs into SOGs (via BFS/Connected Components)
        return consolidate_into_sogs(refined_pairs)


class SOG(Orthogroup):
    def __init__(self, sog_id: str | None = None, genes: list[Gene] | None = None) -> None:
        """Initialize an SOG.

        Args:
            og_id (str): The unique identifier for the sog.
            genes (list[Gene]): The genes in the sog.
        """
        super().__init__(sog_id, genes)

    def add_gene(self, gene: Gene):
        """Add a gene to the SOG.

        If the gene is not already present in the SOG, it will be added and its og attribute set to this SOG.

        Args:
            gene (Gene): The gene to add.
        """
        if gene not in self._genes:
            gene.sog = self
            self._genes.append(gene)

    def refine(self, *args, **kwargs):
        """
        Disable this because a SOG is already refined.
        It shouldn't be refined again.
        """
        raise NotImplementedError("SOG objects are already refined and cannot be refined further.")


def compare_gene_sets(
    genes_a: list[Gene], genes_b: list[Gene], window_size: int, ratio_threshold: float
) -> list[tuple[Gene, Gene]]:
    """Compare two sets of genes from different genomes within the same Orthogroup and return a list of syntenically supported
    orthologous pairs.

    The comparison is based on the given window size and ratio threshold. Genes are considered syntenic if they share enough
    overlapping orthogroups within the specified window, as determined by the ratio threshold.

    Args:
        genes_a (list[Gene]): A list of genes from the first genome in the Orthogroup.
        genes_b (list[Gene]): A list of genes from the second genome in the Orthogroup.
        window_size (int): The size of the window used for synteny comparison. This determines how many genes around a focal gene
            are considered when comparing synteny.
        ratio_threshold (float): The threshold for the ratio of overlapping genes to consider two genes as syntenic. If the ratio
            of overlapping orthogroups between two genes is greater than or equal to this threshold, they are considered syntenic.

    Returns:
        list[tuple[Gene, Gene]]: A list of tuples where each tuple contains a pair of syntenically supported orthologous genes
        from the two input gene sets.
    """
    # Get the shared orthogroups once for this pair
    genome_a: Genome = genes_a[0].genome
    genome_b: Genome = genes_b[0].genome
    shared_ogs = get_shared_ogs(genome_a, genome_b)

    # Determine the primary and secondary genomes
    primary, secondary = (genes_a, genes_b) if len(genes_a) <= len(genes_b) else (genes_b, genes_a)

    # Pre-calculate all windows for the secondary set
    secondary_windows = {
        gene: [g.og.id for g in genome_b.get_window(gene, shared_ogs, window_size) if g.og] for gene in secondary
    }

    refined_pairs = []

    # Compare windows
    for focal_gene in primary:
        genes = focal_gene.genome.get_window(focal_gene, shared_ogs, window_size)
        win_p = [g.og.id for g in genes if g.og]

        # Find the best match
        best_candidate = None
        max_r = -1.0

        for candidate, win_s in secondary_windows.items():
            ratio = calculate_synteny_ratio(win_p, win_s)

            if ratio >= ratio_threshold and ratio > max_r:
                max_r = ratio
                best_candidate = candidate

        if best_candidate:
            refined_pairs.append((focal_gene, best_candidate))

    return refined_pairs


def get_shared_ogs(genome_a: Genome, genome_b: Genome) -> set[str]:
    """Retrieve the set of Orthogroup IDs that are present in both provided genomes.

    This function is used to identify the common orthogroups shared between two genomes, which serves as a basis for comparing
    their synteny and identifying syntenic orthologous pairs.

    Args:
        genome_a (Genome): The first genome.
        genome_b (Genome): The second genome.

    Returns:
        set[str]: A set of Orthogroup IDs that are common to both genomes.
    """
    ogs_a = {g.og.id for g in genome_a if g.og}
    ogs_b = {g.og.id for g in genome_b if g.og}

    return ogs_a.intersection(ogs_b)


def align_sog_dict(sog_dict: dict[Gene, list[Gene]]) -> dict[Gene, list[Gene | None]]:
    """Aligns the neighborhood lists within the dictionary so that the focal gene (the key) is at the same index in every value.

    The function pads the neighborhoods with 'None' values to ensure that the focal gene is always at the same index across all
    neighborhoods. This alignment is necessary for further analysis, such as synteny comparison.

    Args:
        sog_dict (dict[Gene, list[Gene]]): A dictionary where keys are genes and values are lists of neighboring genes.

    Returns:
        dict[Gene, list[Gene | None]]: A dictionary with the same structure as `sog_dict`, but with padded neighborhoods to align
        the focal gene at the same index.
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


def calculate_synteny_ratio(win_a: list[str], win_b: list[str]) -> float:
    """Calculates the 1-to-1 synteny match ratio between two dynamic windows.

    Specifically, this function computes the ratio of overlapping orthogroups present in both window sets. The overlap is
    determined by finding the minimum count for each shared Orthogroup ID across the two windows.

    Args:
        win_a (list[str]): A list of Orthogroup IDs representing the first dynamic window.
        win_b (list[str]): A list of Orthogroup IDs representing the second dynamic window.

    Returns:
        float: The synteny ratio, calculated as the number of overlapping orthogroups divided by the length of the longer window.
        If either window is empty, the function returns 0.0.
    """
    if not win_a or not win_b:
        return 0.0

    # Counter handles the 1-to-1 matching via the & (intersection) operator
    counts_a = Counter(win_a)
    counts_b = Counter(win_b)

    # This automatically takes the minimum count for each shared ID
    matches = sum((counts_a & counts_b).values())

    return matches / max(len(win_a), len(win_b))


def consolidate_into_sogs(pairs: list[tuple[Gene, Gene]]) -> list[SOG]:
    """Consolidates a list of gene pairs into Syntenic Orthologous Groups (SOGs).

    This function uses the provided list of syntenically supported orthologous pairs to construct connected components within a
    graph representation. Each connected component represents a SOG, where all genes in the component are syntenically linked.

    Args:
        pairs (list[tuple[Gene, Gene]]): A list of tuples, each containing two syntenically supported orthologous genes from
        different genomes.

    Returns:
        list[SOG]: A list of SOG (refined orthogroup) objects.
    """
    # Build an adjacency list (the graph)
    adj = defaultdict(list)
    nodes = set()
    for u, v in pairs:
        adj[u].append(v)
        adj[v].append(u)
        nodes.update([u, v])

    visited = set()
    sog_clusters: list[SOG] = []

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
            sog_clusters.append(SOG(genes=component))

    return sog_clusters


def _perform_pairwise_comparisons(orthogroup: Orthogroup, window_size: int, ratio_threshold: float) -> list[tuple[Gene, Gene]]:
    """A helper function to organize genes by genome and run comparisons.

    Args:
        window_size (int): The size of the window used for synteny comparison.
        ratio_threshold (float): The threshold for the ratio of overlapping genes to consider two genes as syntenic.

    Returns:
        list[tuple[Gene, Gene]]: A list of refined pairs of genes that are considered syntenic based on the given criteria.
    """
    genes_by_genome = defaultdict(list)
    for gene in orthogroup:
        genes_by_genome[gene.genome].append(gene)

    all_refined_pairs = []
    # Combinations ensures we only compare Genome A to Genome B once
    for genome_a, genome_b in combinations(genes_by_genome.keys(), 2):
        pairs = compare_gene_sets(genes_by_genome[genome_a], genes_by_genome[genome_b], window_size, ratio_threshold)
        all_refined_pairs.extend(pairs)
    return all_refined_pairs
