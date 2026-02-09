from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ._orthogroup import Orthogroup


class Protein:
    def __init__(self, protein_id: str):
        self.protein_id = protein_id

        # For pointers
        self.gene: Gene | None = None


class Proteome:
    def __init__(self, sample_name: str):
        self.sample_name = sample_name
        self._proteins: dict[str, Gene] = {}

    def __len__(self) -> int:
        """Returns the number of genes in the genome."""
        return len(self._proteins)

    def __getitem__(self, key: str) -> Gene:
        """Returns the gene at the given index."""
        return self._proteins[key]

    def add_protein(self, protein_obj: Protein):
        """Adds a gene object to the genome.
        Args:
            gene_obj (Gene): The Gene object to add.
        """
        self._proteins[protein_obj.protein_id] = protein_obj


class Gene:
    """Represents a gene in a genome.

    Attributes:
        locus_tag (str): The locus tag of the gene.
        prod_acc (str): The product accession number of the gene.
        scaffold (str): The scaffold on which the gene is located.
        operon_id (int): The ID of the operon the gene belongs to (default: -1).
        genome (Genome | None): The Genome object this gene belongs to.
        orthogroup (Orthogroup | None): The Orthogroup object representing the gene's orthogroup.
        chain_index (int | None): The chain index of the gene.
    """

    __slots__ = ("seqid", "id", "start", "end", "proteins", "operon_id", "genome", "orthogroup", "chain_index", "len")

    def __init__(self, seqid: str, start: int, end: int, id: str):
        self.seqid = seqid
        self.start = start
        self.len = abs(end - start)
        self.id = id

        # For operon detection (if needed)
        self.operon_id = -1

        # For pointers
        self.genome: Genome | None = None
        self.proteins: list[Protein] = []
        self.orthogroup: Orthogroup | None = None
        self.chain_index: int | None = None

    def __repr__(self):
        og_id = self.orthogroup.og_id if self.orthogroup else "None"
        return f"[{self.id} | {og_id}]"


class Genome:
    """Represents a genome.

    Attributes:
        sample_name (str): The name of the sample.
        chromosome_type (str): The type of chromosome ('c' for circular, 'l' for linear).
        _genes (list[Gene]): A list of Gene objects in the genome.
    """

    __slots__ = ("sample_name", "chromosome_type", "_genes", "_gene_map")

    def __init__(self, sample_name: str, chromosome_type: str = "l"):
        self.sample_name = sample_name
        self.chromosome_type = chromosome_type
        self._genes: list[Gene] = []
        self._gene_map: dict[str, Gene] = {}

    def __len__(self) -> int:
        """Returns the number of genes in the genome."""
        return len(self._genes)

    def __getitem__(self, key: int | str) -> Gene:
        """Returns the gene at the given index."""
        if isinstance(key, str):
            return self._gene_map[key]
        return self._genes[key]

    def __iter__(self):
        """Returns an iterator for the genes in the genome."""
        return iter(self._genes)

    def add_gene(self, gene_obj: Gene):
        """Adds a gene object to the genome.
        Args:
            gene_obj (Gene): The Gene object to add.
        """
        gene_obj.genome = self
        gene_obj.chain_index = len(self._genes)
        self._genes.append(gene_obj)
        self._gene_map[gene_obj.id] = gene_obj

    def get_window(self, focal_gene: Gene, target_ogs: set, window_size: int) -> list[str]:
        """
        Builds a neighborhood using ONLY OG IDs found in the target_ogs set.
        Args:
            focal_gene: The Gene to center the window around.
            target_ogs: A set of Orthogroup IDs to include in the neighborhood.
            window_size: The size of the window to build around the focal gene.

        Returns:
            A list of Orthogroup IDs found in the target_ogs set within the window.
        """
        half_win = window_size // 2
        neighborhood = []
        total_genes = len(self)

        for direction in [-1, 1]:
            found_count, offset = 0, 1
            while found_count < half_win:
                curr_idx = focal_gene.chain_index + (offset * direction)

                # Boundary checks
                if self.chromosome_type == "l" and (curr_idx < 0 or curr_idx >= total_genes):
                    break

                neighbor_gene = self[curr_idx % total_genes]
                if neighbor_gene.seqid != focal_gene.seqid:
                    break

                # Logic is now: Is this gene a shared anchor?
                if neighbor_gene.orthogroup and neighbor_gene.orthogroup.og_id in target_ogs:
                    neighborhood.append(neighbor_gene.orthogroup.og_id)
                    found_count += 1

                offset += 1
                if offset >= total_genes:
                    break
        return neighborhood
