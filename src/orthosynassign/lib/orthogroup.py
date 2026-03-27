from __future__ import annotations

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
