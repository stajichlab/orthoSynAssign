from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ._orthogroup import Orthogroup


class Gene:
    """
    Represents a gene in a genome.
    """

    __slots__ = ("locus_tag", "prod_acc", "scaffold", "operon_id", "genome", "orthogroup", "chain_index")

    def __init__(self, locus_tag: str, prod_acc: str, scaffold: str):
        self.locus_tag = locus_tag
        self.prod_acc = prod_acc
        self.scaffold = scaffold

        # For operon detection (if needed)
        self.operon_id = -1

        # For triple pointers
        self.genome: Genome | None = None
        self.orthogroup: Orthogroup | None = None
        self.chain_index: int | None = None

    def __repr__(self):
        og_id = self.orthogroup.og_id if self.orthogroup else "None"
        return f"[{self.locus_tag} | {og_id}]"


class Genome:
    __slots__ = ("sample_name", "chromosome_type", "_genes")

    def __init__(self, sample_name: str, chromosome_type: str = "c"):
        self.sample_name = sample_name
        self.chromosome_type = chromosome_type
        self._genes: list[Gene] = []

    def __len__(self) -> int:
        return len(self._genes)

    def __getitem__(self, index: int) -> Gene:
        return self._genes[index]

    def __iter__(self):
        return iter(self._genes)

    def add_gene(self, gene_obj: Gene):
        gene_obj.genome = self
        gene_obj.chain_index = len(self._genes)
        self._genes.append(gene_obj)

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
                if neighbor_gene.scaffold != focal_gene.scaffold:
                    break

                # Logic is now: Is this gene a shared anchor?
                if neighbor_gene.orthogroup and neighbor_gene.orthogroup.og_id in target_ogs:
                    neighborhood.append(neighbor_gene.orthogroup.og_id)
                    found_count += 1

                offset += 1
                if offset >= total_genes:
                    break
        return neighborhood
