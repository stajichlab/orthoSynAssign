# rs.pyi
from __future__ import annotations

class SyntenyEngine:
    """
    High-performance Rust backend for synteny analysis.
    """
    def __init__(self, num_ogs: int, ogs_all: list[list[int]], seqids_all: list[list[int]], is_circular_all: list[bool]) -> None:
        """Initializes a SyntenyEngine with genomic data.

        Args:
            num_ogs (int): The total number of orthogroups.
            ogs_all (list[list[int]]): A list of lists of orthogroup indices for genes in each genome.
            seqids_all (list[list[int]]): A list of lists of sequence/scaffold indices for genes in each genome.
            is_circular_all (list[bool]): A list of booleans indicating whether each genome is circular or not.
        """
        ...

    def refine(self, og_idx: int, window_size: int, ratio_threshold: float) -> list[list[tuple[int, int]]]:
        """Coordinates the refinement of a single Orthogroup using physical anchors.

        Args:
            og_idx (int): The index of the orthogroup to refine.
            window_size (int): The size of the window to build around the genes in the orthogroup.
            ratio_threshold (float): The minimum ratio to consider for synteny.

        Returns:
            list[list[tuple[int, int]]]: A list of clusters, where each cluster is a list of (genome_idx, gene_idx) physical
                anchors.
        """
        ...

class VisualizeEngine:
    """
    High-performance Rust backend for visualization.
    """
    def __init__(
        self,
        sogs: list[list[tuple[int, int]]],
        ogs_all: list[list[int]],
        seqids_all: list[list[int]],
        is_circular_all: list[bool],
    ) -> None:
        """Initializes a VisualizeEngine with genomic data.

        Args:
            sogs (list[list[tuple[int, int]]]): A list of gene indices (genome_idx, gene_idx) for all refined orthogroups.
            ogs_all (list[list[int]]): A list of lists of orthogroup indices for genes in each genome.
            seqids_all (list[list[int]]): A list of lists of sequence/scaffold indices for genes in each genome.
            is_circular_all (list[bool]): A list of booleans indicating whether each genome is circular or not.
        """
        ...

    def get_aligned_og(
        self, sog_idx: int, window_size: int, keep_all_genes: bool = False
    ) -> list[tuple[tuple[int, int], list[int | None]]]:
        """Retrieves a list of genes and their corresponding windows aligned by the genes from the given orthogroup.

        Args:
            sog_idx (int): The index of the refined orthogroup to visualize.
            window_size (int): The size of the window to build around the genes in the orthogroup.
            keep_all_genes (bool): whether to keep all genes even without orthogroup assignment.

        Returns:
            list[tuple[tuple[int, int], list[int | None]]]: A list of tuples where the first element is the focal gene
                indices and the second is the aligned window of gene indices (None for padding).
        """
        ...

def get_window(
    seqid_vec: list[int], og_vec: list[int], gene_idx: int, target_ogs: list[int], window_size: int, is_circular: bool
) -> list[int]:
    """Retrieve the neighborhood gene indices from the focal gene index with a given window size.

    Args:
        seqid_vec (list[int]): An array of sequence/scaffold IDs for genes in the genome.
        og_vec (list[int]): An array of orthogroup indices for genes in the genome.
        gene_idx (int): The index of the focal gene.
        target_ogs (list[int]): A list of target orthogroups to consider when building the window.
        window_size (int): The size of the window to build around the focal gene.
        is_circular (bool): A boolean indicating whether the genome is circular or not.

    Returns:
        list[int]: A list of gene indices found within the window.
    """
    ...

def calculate_synteny_ratio(win_a: list[int], win_b: list[int]) -> float:
    """Calculates the 1-to-1 synteny match ratio between two dynamic windows.

    Args:
        win_a (list[int]): A list of orthogroup indices representing the first dynamic window.
        win_b (list[int]): A list of orthogroup indices representing the second dynamic window.

    Returns:
        float: The synteny ratio, calculated as the number of overlapping orthogroups divided by the length of the
            longer window. If either window is empty, returns 0.0.
    """
    ...
