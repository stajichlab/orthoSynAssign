from __future__ import annotations

from typing import TYPE_CHECKING, Any, Iterator, Literal

if TYPE_CHECKING:
    from .orthogroup import Orthogroup


class Gene:
    """Represents a gene in a genome.

    Attributes:
        seqid (str): The chromosome/scaffold/contig on which the gene is located.
        start (int): The start position of the gene.
        len (int): The length of the gene.
        id (str): The name of the gene.
        genome (Genome | None): The genome to which the gene belongs, if any.
        orthogroup (Orthogroup | None): The orthogroup to which the gene belongs, if any.
        index (int | None): The index of the gene in its genome.
    """

    __slots__ = ("seqid", "start", "len", "id", "genome", "orthogroup", "index")

    def __init__(self, seqid: str, start: int, end: int, gene_id: str) -> None:
        """Initialize a new Gene object.

        Args:
            seqid (str): The chromosome/scaffold/contig on which the gene is located.
            start (int): The start position of the gene.
            end (int): The end position of the gene.
            gene_id (str): The name of the gene.
        """
        self.seqid = seqid
        self.start = start
        self.len = abs(end - start)
        self.id = gene_id

        # For operon detection (if needed)
        # self.operon_id = -1

        # For pointers
        self.genome: Genome | None = None
        self.orthogroup: Orthogroup | None = None
        self.index: int | None = None

    def __repr__(self) -> str:
        """Return a string representation of the Gene object.

        Returns:
            str: A string containing the gene name, genome name, and orthogroup ID.
        """
        genome_name = self.genome.name if self.genome else "Unknown genome"

        og_id = self.orthogroup.id if self.orthogroup else "Unassigned orthogroup"

        return f"{self.id} @ {genome_name} | {og_id}"

    def __getstate__(self) -> dict[str, Any]:
        """Return the state of the Gene object as a dictionary.

        Returns:
            dict[str, Any]: A dictionary containing all attributes of the Gene object.
        """
        return {slot: getattr(self, slot) for slot in self.__slots__ if hasattr(self, slot)}

    def __setstate__(self, state: dict[str, Any]) -> None:
        """Set the state of the Gene object from a dictionary.

        Args:
            state (dict[str, Any]): A dictionary containing all attributes of the Gene object.
        """
        for slot, value in state.items():
            setattr(self, slot, value)


class Genome:
    """Represents a genome.

    Attributes:
        name (str): The name of the sample.
        chromosome_type (str): The type of chromosome ('c' for circular, 'l' for linear).
        _genes (list[Gene]): A list to store all gene objects associated with this genome.
        _gene_map (dict[str, Gene]): A dictionary that maps gene IDs to their corresponding Gene objects for fast lookup.
    """

    __slots__ = ("name", "chromosome_type", "_genes", "_gene_map")

    def __init__(self, name: str, *, chromosome_type: Literal["l", "c"] = "l") -> None:
        """Initialize a new Genome object.

        Args:
            name (str): The name of the genome.
            chromosome_type (Literal["l", "c"], optional): The type of chromosome ('c' for circular, 'l' for linear). Defaults to
            'l'.

        Raises:
            ValueError: If an invalid chromosome type is provided.
        """
        self.name = name
        if chromosome_type not in ("l", "c"):
            raise ValueError("Invalid chromosome type. Must be 'l' or 'c'.")
        self.chromosome_type = chromosome_type
        self._genes: list[Gene] = []
        self._gene_map: dict[str, Gene] = {}

    def __repr__(self) -> str:
        """Return a string representation of the Genome object.

        Returns:
            str: A string containing the genome name and type.
        """
        chromosome_type = "linear" if self.chromosome_type == "l" else "circular"
        return f"[{self.name} ({chromosome_type}) | with {len(self._genes)} genes]"

    def __len__(self) -> int:
        """Returns the number of genes in the genome.

        Returns:
            int: The total number of genes in the genome.
        """
        return len(self._genes)

    def __getitem__(self, key: int | str | slice) -> Gene | list[Gene]:
        """Returns the gene at the given index or by ID.

        Args:
            key (int | str | slice): The index, ID, or slice of genes to retrieve.

        Returns:
            Gene | list[Gene]: A Gene object if a single item is requested, otherwise a list of Gene objects.

        Raises:
            KeyError: If the provided string ID does not correspond to any gene in the genome.
        """
        if isinstance(key, str):
            return self._gene_map[key]
        return self._genes[key]

    def __iter__(self) -> Iterator[Gene]:
        """Returns an iterator for the genes in the genome.

        Returns:
            Iterator[Gene]: An iterator object that allows iteration over all genes in the genome.
        """
        return iter(self._genes)

    def __getstate__(self) -> dict[str, Any]:
        """Return the state of the Genome object as a dictionary.

        Returns:
            dict[str, Any]: A dictionary containing all attributes of the Genome object.
        """
        return {slot: getattr(self, slot) for slot in self.__slots__ if hasattr(self, slot)}

    def __setstate__(self, state: dict[str, Any]) -> None:
        """Set the state of the Genome object from a dictionary.

        This method ensures that back-pointers to the genome are correctly restored when deserializing the object.

        Args:
            state (dict[str, Any]): A dictionary containing all attributes of the Genome object.
        """
        for slot, value in state.items():
            setattr(self, slot, value)

        # Double check: ensure back-pointers are still correct
        for gene in self._genes:
            gene.genome = self

    def add_gene(self, gene_obj: Gene) -> None:
        """Adds a gene object to the genome.

        This method updates the genome attribute of the gene object and adds it to both internal lists.

        Args:
            gene_obj (Gene): The Gene object to add.
        """
        gene_obj.genome = self
        gene_obj.index = len(self._genes)
        self._genes.append(gene_obj)
        self._gene_map[gene_obj.id] = gene_obj

    def get_window(self, focal_gene: Gene, target_ogs: set[str], window_size: int) -> list[Gene]:
        """Builds a neighborhood using ONLY OG IDs found in the target_ogs set.

        This method searches for genes in both upstream and downstream directions, respecting chromosome type boundaries.

        Args:
            focal_gene (Gene): The Gene to center the window around.
            target_ogs (set[str]): A set of Orthogroup IDs to include in the neighborhood.
            window_size (int): The size of the window to build around the focal gene.

        Returns:
            list[Gene]: A list of ordered Gene objects found in the target_ogs set within the window.
        """
        half_win: int = window_size // 2
        upstream: list[Gene] = []
        downstream: list[Gene] = []
        total_genes: int = len(self._genes)

        for direction in [-1, 1]:
            found_count, offset = 0, 1
            while found_count < half_win:
                curr_idx = focal_gene.index + (offset * direction)

                # Boundary checks for linear chromosomes
                if self.chromosome_type == "l" and (curr_idx < 0 or curr_idx >= total_genes):
                    break

                # Handle circularity with modulo
                neighbor_gene = self._genes[curr_idx % total_genes]

                # Stop if we hit a different scaffold/chromosome
                if neighbor_gene.seqid != focal_gene.seqid:
                    break

                # Collect the gene
                if direction == -1:
                    upstream.append(neighbor_gene)
                else:
                    downstream.append(neighbor_gene)

                # Check if this gene is an "anchor" (exists in target_ogs)
                if neighbor_gene.orthogroup and neighbor_gene.orthogroup.id in target_ogs:
                    found_count += 1

                offset += 1
                if offset >= total_genes:
                    break

        # We combine upstream (reversed) and downstream, excluding the focal gene
        flanks: list[Gene] = upstream[::-1] + downstream
        return [g for g in flanks if g.orthogroup and g.orthogroup.id in target_ogs]
