# rs.pyi
from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import numpy as np
    from numpy.typing import NDArray

class SyntenyEngine:
    """
    High-performance Rust backend for synteny analysis.
    """
    def __init__(
        self, num_orthogroups: int, og_list_all: list[list[int]], seqid_list_all: list[list[int]], is_circular_all: list[bool]
    ) -> None:
        """Initializes a SyntenyEngine with genomic data.

        Args:
            num_orthogroups (int): The total number of unique orthogroups.
            og_inputs (list[NDArray[np.int16]]): A list of arrays of orthogroup IDs for genes in each genome.
            seq_inputs (list[NDArray[np.int32]]): A list of arrays of sequence/scaffold IDs for genes in each genome.
        """
        ...

    def refine(self, og_idx: int, window_size: int, ratio_threshold: float) -> list[list[tuple[int, int]]]:
        """Coordinates the refinement of a single Orthogroup using physical anchors.

        Returns:
            list[list[tuple[int, int]]]: A list of clusters, where each cluster is a list of (genome_idx, gene_idx) physical
                anchors.
        """
        ...

def get_window(og_masked_array: NDArray[np.bool_], seq_array: NDArray[np.int16], gene_idx: int, window_size: int) -> list[int]:
    """Retrieve the neighborhood gene indices from the focal gene index with a given window size.

    This function is used to find the genes within a specified window around a focal gene. It takes into account the orthogroup
    mask and sequence array to identify relevant genes.

    Args:
        og_masked_array (NDArray[np.bool_]): A boolean array representing whether each gene belongs to the target orthogroup.
        seq_array (NDArray[np.int16]): An array of sequence/scaffold IDs for genes in the genome.
        gene_idx (int): The index of the focal gene.
        window_size (int): The size of the window to build around the focal gene.

    Returns:
        list[int]: A list of gene indices found within the window.
    """
    ...

def calculate_synteny_ratio(win_a: list[int], win_b: list[int]) -> float:
    """Calculates the 1-to-1 synteny match ratio between two dynamic windows.

    Specifically, this function computes the ratio of overlapping orthogroups present in both window sets. The overlap is
    determined by finding the minimum count for each shared Orthogroup ID across the two windows.

    Args:
        win_a (list[int]): A list of Orthogroup indices representing the first dynamic window.
        win_b (list[int]): A list of Orthogroup indices representing the second dynamic window.

    Returns:
        float: The synteny ratio, calculated as the number of overlapping orthogroups divided by the length of the longer window.
        If either window is empty, the function returns 0.0.
    """
    ...
