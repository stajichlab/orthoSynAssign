"""
File parsers for GTF, GFF3, and OrthoFinder ortholog tables.
"""

from __future__ import annotations

import gzip
import logging
import re
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any

from ._gene import Gene, Genome, Protein, Proteome
from ._orthogroup import Orthogroup

logger = logging.getLogger(__name__)


class AnnotationParser(ABC):
    """Abstract base class for genomic annotation parsers."""

    def __init__(self, file: str | Path):
        self.file = Path(file)
        if not self.file.exists():
            raise FileNotFoundError(f"Annotation file not found: {self.file}")
        if not self.file.is_file():
            raise ValueError(f"Annotation file must be a regular file: {self.file}")
        self._genome = None
        self._proteome = None

    def parse(self, *args: Any, **kwargs: Any):
        """Parse annotation file."""
        try:
            opener = gzip.open if self._is_gzip_file() else open

            with opener(self.file, "rt", encoding="utf-8") as f:
                logger.info(f"Reading annotation file: {self.file}")
                sample = re.sub(r"\.(bed|gff3?|gtf)(\.gz)?$", "", self.file.name)
                self._genome = Genome(sample_name=sample)
                self._proteome = Proteome(sample_name=sample)

                for line_num, line in enumerate(f, 1):
                    line = line.strip()

                    # Skip empty lines and comments
                    if not line or line.startswith("#"):
                        continue

                    self._parser(line_num, line, *args, **kwargs)

            logger.info(f"Loaded {len(self._genome)} features from {self.file}")
        except Exception:
            self._genome = None
            self._proteome = None
            raise

        return self._genome

    @abstractmethod
    def _parser(self, line_num: int, line: str, *args: Any, **kwargs: Any):
        pass

    def _is_gzip_file(self) -> bool:
        """Checks whether a file is in gzip format.

        This function reads the first two bytes of the file and checks if they match the gzip magic number.

        Args:
            file (str | Path): The path to the file to check.

        Returns:
            bool: True if the file is in gzip format, False otherwise.
        """
        with open(self.file, "rb") as f:
            magic_number = f.read(2)
        return magic_number == b"\x1f\x8b"


class BedParser(AnnotationParser):
    def _parser(self, line_num: int, line: str):
        fields = line.split("\t")
        if len(fields) != 4:
            raise ValueError(f"Line {line_num}: Expected 4 fields, got {len(fields)}")
        (seqid, start, end, name) = fields
        self._genome.add_gene(Gene(seqid=seqid, start=int(start), end=int(end), id=name))


class _GFFGTFParser(AnnotationParser):
    """
    Abstract base class for GFF and GTF parsers.
    """

    def parse(
        self,
        gene_parser: dict[str, str] = {"feature_type": ["gene"], "ID": "locus_tag"},
        protein_parser: dict[str, str] = {"feature_type": ["CDS"], "ID": "protein_id"},
        *args,
        **kwargs,
    ):
        self._g_par = gene_parser
        self._p_par = protein_parser
        self._last_gene = None
        return super().parse(*args, **kwargs)

    def _parser(
        self,
        line_num: int,
        line: str,
        *args,
        **kwargs,
    ):
        fields = line.split("\t")
        if len(fields) != 9:
            raise ValueError(f"Line {line_num}: Expected 9 fields, got {len(fields)}")
        (
            seqid,
            _,
            feature_type,
            start,
            end,
            _,
            _,
            _,
            attributes,
        ) = fields

        try:
            if feature_type in self._g_par["feature_type"]:
                id = self._attrs_parser(attributes, parsed_key=self._g_par["ID"])
                gene = Gene(seqid=seqid, id=id, start=start, end=end)
                self._genome.add_gene(gene)
                self._last_gene = gene
            elif feature_type in self._p_par["feature_type"]:
                id = self._attrs_parser(attributes, parsed_key=self._p_par["ID"])
                protein = Protein(protein_id=id)
                protein.gene = self._last_gene if self._last_gene else None
                self._proteome.add_protein(protein)

        except ValueError as e:
            raise ValueError(f"Line {line_num}: Error parsing attributes - {e}")

    @abstractmethod
    def _attrs_parser(self, attr_string, parsed_key):
        """
        Abstract method to parse the attributes string for a specific file format.
        Subclasses must implement this method to parse the attributes.
        """
        pass


class GFFParser(_GFFGTFParser):
    """
    Parse GFF3 attributes field.

    GFF3 format: ID=gene1;Name=MyGene;Note=Some note

    Args:
        attr_string: The attributes column from a GFF3 line

    Returns:
        Dictionary of attribute key-value pairs
    """

    def _attrs_parser(self, attr_string, parsed_key):
        if not attr_string or attr_string == ".":
            return None

        for item in attr_string.split(";"):
            item = item.strip()
            if "=" in item:
                key, value = item.split("=", 1)
                if key == parsed_key:
                    return value
        return None


class GTFParser(_GFFGTFParser):
    """
    Parse GTF attributes field.

    GTF format: gene_id "ENSG00000123"; transcript_id "ENST00000456";

    Args:
        attr_string: The attributes column from a GTF line

    Returns:
        Dictionary of attribute key-value pairs
    """

    def _attrs_parser(self, attr_string, parsed_key):
        # Split by semicolon and process each key-value pair
        for item in attr_string.split(";"):
            item = item.strip()
            if item:
                parts = item.split(None, 1)  # Split on first whitespace
                if len(parts) == 2:
                    key = parts[0]
                    value = parts[1].strip('"')  # Remove quotes
                    if key == parsed_key:
                        return value
        return None


def read_orthofinder_table(file: str, genomes: dict[str, Genome]) -> list[Orthogroup]:
    """
    Read an OrthoFinder Orthogroups.tsv file.

    The file format is tab-separated with:
    - First column: Orthogroup ID
    - Subsequent columns: Gene IDs for each species (column headers are species names)

    Args:
        file_path: Path to the Orthogroups.tsv file

    Returns:
        Dictionary mapping orthogroup IDs to dictionaries of species->gene lists
        Example:
        {
            'OG0000001': {
                'Species1': ['gene1', 'gene2'],
                'Species2': ['gene3']
            }
        }

    Raises:
        FileNotFoundError: If the file doesn't exist
        ValueError: If the file format is invalid
    """
    file: Path = Path(file)

    orthogroups = []
    logger.info(f"Reading Orthogroup file: {file}")

    with open(file, "r", encoding="utf-8") as f:
        # Read header to get species names
        header = f.readline().strip().split("\t")
        if len(header) < 3:
            raise ValueError("Invalid Orthogroup file: expected at least 3 columns")

        # First column is "Orthogroup", rest are sample names
        samples = header[1:]
        logger.info("Found %s samples: %s", len(samples), ", ".join(samples))
        if genomes:
            sample_set = set(header[1:])
            genome_keys = set(genomes.keys())

            # Find the difference
            missing = sample_set - genome_keys

            if missing:
                raise ValueError(
                    f"The following samples were not found in loaded genomes: {', '.join(missing)}\n"
                    "Please check whether the samples in the orthogroup file match to the gff files.\n"
                    f"{genome_keys}"
                )

        # Read orthogroup data
        for line_num, line in enumerate(f, 2):  # Start at 2 since we read header
            line = line.strip("\n")
            if not line:
                continue

            fields = line.split("\t")
            if len(fields) != len(header):
                logger.warning(f"Line {line_num}: Expected {len(header)} fields, got {len(fields)}")
                continue

            orthogroup = Orthogroup(fields[0])

            # Parse genes for each sample
            for i, sample in enumerate(samples, 1):
                if genomes:
                    genome = genomes[sample]
                else:
                    genome = Genome(sample)
                genes_str = fields[i].strip()
                if genes_str:
                    # Genes are usually separated by commas or spaces
                    genes = [g.strip() for g in genes_str.replace(",", " ").split() if g.strip()]
                    for gene in genes:
                        orthogroup.add_gene(genome[gene])

            orthogroups.append(orthogroup)

    logger.info(f"Loaded {len(orthogroups)} orthogroups from {file}")
    return orthogroups


def save_results_tsv(results: list[tuple[str, dict[Genome, list[Gene]]]], filename):
    """
    Format: SOG_ID \t Genome_A \t Genome_B \t ...
    Where cells contain comma-separated locus_tags.
    """
    if not results:
        return

    # Determine all unique genomes across all results for header
    all_genomes = sorted(list({g.sample_name for _, sog in results for g in sog.keys()}))

    with open(filename, "w") as f:
        # Write Header
        header = ["SOG_ID"] + all_genomes
        f.write("\t".join(header) + "\n")

        # Write Rows
        for sog_id, sog_dict in results:
            row = [sog_id]
            for g_name in all_genomes:
                # Find the genome object that matches the name
                # (Assuming sog_dict keys are Genome objects)
                matching_genes = []
                for genome_obj, genes in sog_dict.items():
                    if genome_obj.sample_name == g_name:
                        matching_genes = [g.id for g in genes]
                        break

                row.append(",".join(matching_genes) if matching_genes else "")
            f.write("\t".join(row) + "\n")
