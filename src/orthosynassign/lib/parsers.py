"""
File parsers for GTF, GFF3, and OrthoFinder ortholog tables.
"""

from __future__ import annotations

import gzip
import logging
import re
from abc import ABC, abstractmethod
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterator

from .gene import Gene, Genome
from .orthogroup import SOG, Orthogroup

logger = logging.getLogger(__name__)


class AnnotationParser(ABC):
    """Abstract base class for genomic annotation parsers.

    This class provides a framework for parsing various types of genomic annotation files such as GTF,
    GFF3, and others. It supports reading compressed files in gzip format.

    Attributes:
        file (Path): The path to the annotation file.
        _genome (Genome | None): A `Genome` object representing the parsed data from the file.
            This attribute is initialized to `None` and populated during the parsing process.
    """

    __slots__ = ("file", "_genome")

    def __init__(self, file: str | Path) -> None:
        """Initialize an AnnotationParser instance.

        Args:
            file (str | Path): The path to the annotation file.

        Raises:
            FileNotFoundError: If the specified annotation file does not exist.
            ValueError: If the specified annotation file is not a regular file.
        """
        self.file = Path(file)
        if not self.file.exists():
            raise FileNotFoundError(f"Annotation file not found: {self.file}")
        if not self.file.is_file():
            raise ValueError(f"Annotation file must be a regular file: {self.file}")
        self._genome: Genome | None = None

    def parse(self, *args: Any, **kwargs: Any) -> Genome:
        """Parse annotation file.

        This method reads the annotation file line by line, skipping empty lines and comments. It delegates the actual parsing of
        each line to the `_parser` method, which must be implemented by subclasses. After processing all lines, it logs the number
        of features loaded into the `Genome` object.

        Args:
            *args (Any): Additional positional arguments to pass to the `_parser` method.
            **kwargs (Any): Additional keyword arguments to pass to the `_parser` method.

        Returns:
            Genome: The `Genome` object containing the parsed data from the file.

        Raises:
            FileNotFoundError: If the specified annotation file does not exist.
            ValueError: If the specified annotation file is not a regular file.
            Exception: Any exception raised during the parsing process, which will result in resetting the `_genome` attribute to
                `None`.
        """
        try:
            opener = gzip.open if self._is_gzip_file() else open

            with opener(self.file, "rt", encoding="utf-8") as f:
                logger.info(f"Reading annotation file: {self.file}")
                sample = re.sub(r"\.(bed|gff3?|gtf)(\.gz)?$", "", self.file.name)
                self._genome = Genome(sample)

                for line_num, line in enumerate(f, 1):
                    line = line.strip()

                    # Skip empty lines and comments
                    if not line or line.startswith("#"):
                        continue

                    self._parser(line_num, line, *args, **kwargs)

            logger.info(f"Loaded {len(self._genome)} features from {self.file}")
        except Exception:
            self._genome = None
            raise

        return self._genome

    @abstractmethod
    def _parser(self, line_num: int, line: str, *args: Any, **kwargs: Any) -> None:
        """Abstract method for parsing a single line of an annotation file.

        This method must be implemented by subclasses to parse the specific format of the annotation
        file being processed.

        Args:
            line_num (int): The current line number in the annotation file.
            line (str): The content of the current line, stripped of leading and trailing whitespace.
            *args (Any): Additional positional arguments for custom parsing logic.
            **kwargs (Any): Additional keyword arguments for custom parsing logic.

        Raises:
            Exception: Any exception raised during the parsing process will be propagated up to the `parse` method.
        """
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
    """BED parser class for genomic annotation files.

    This class implements the `_parser` method to parse BED format annotation files. Each line in a BED file contains information
    about a genomic feature, such as its chromosomal location and name.
    """

    def _parser(self, line_num: int, line: str) -> None:
        """Parse a single line of a BED format annotation file.

        This method splits the line into fields and extracts the necessary information about the genomic feature,
        such as its chromosomal location and name. It then adds a `Gene` object representing this feature to
        the `Genome` object being built.

        Args:
            line_num (int): The current line number in the annotation file.
            line (str): The content of the current line, stripped of leading and trailing whitespace.

        Raises:
            ValueError: If the line does not contain exactly 4 fields, indicating an invalid format for a BED
                file.
        """
        self._genome: Genome
        fields = line.split("\t")
        if len(fields) != 4:
            raise ValueError(f"Line {line_num}: Expected 4 fields, got {len(fields)}")
        (seqid, start, end, name) = fields
        start, end = int(start), int(end)
        gene = Gene(seqid=seqid, start=start, end=end, gene_id=name)
        self._genome.add_gene(gene)

        names = name.split(";")
        if len(names) > 1:
            for n in names:
                g = Gene(seqid=seqid, start=start, end=end, gene_id=n)
                g.representative = gene
                self._genome.add_gene(g, is_isoform=True)


def read_og_table(file: str | Path, genomes: dict[str, Genome]) -> list[Orthogroup]:
    """Read an OrthoFinder-style orthogroups.tsv file.

    The file format is tab-separated with:
    - First column: Orthogroup ID
    - Subsequent columns: Gene IDs for each species (column headers are species names)

    Args:
        file (str | Path): Path to the orthogroups.tsv file
        genomes (dict[str, Genome]): A dictionary of `Genome` objects keyed by sample name.

    Returns:
        list[Orthogroup]: A list of parsed Orthogroups.
    """
    file: Path = Path(file)
    orthogroups: list[Orthogroup] = []
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
            missing = set(samples) - set(genomes.keys())
            if missing:
                raise ValueError(f"Samples in OG file not found in loaded genomes: {', '.join(missing)}")

        # Read orthogroup data
        for line_num, line in enumerate(f, 2):  # Start at 2 since we read header
            fields = line.strip("\n").split("\t")
            if not fields[0]:
                continue

            if len(fields) != len(header):
                logger.warning(f"Line {line_num}: Expected {len(header)} fields, got {len(fields)}")
                continue

            orthogroup = Orthogroup(fields[0])

            # Parse genes for each sample
            for i, sample in enumerate(samples, 1):
                genome = genomes[sample]
                gene_data = fields[i].strip()
                if not gene_data:
                    continue

                # Genes are usually separated by commas or spaces
                genes = [g.strip() for g in gene_data.replace(",", " ").split()]
                for gene in genes:
                    gene_obj = genome._gene_map.get(gene)
                    if gene_obj:
                        orthogroup.add_gene(gene_obj)
                    else:
                        logger.warning(f"Gene {gene} not found in {sample}")

            orthogroups.append(orthogroup)

    logger.info(f"Loaded {len(orthogroups)} orthogroups from {file}")
    return orthogroups


def write_og_table(results_gen: Iterator[SOG], all_genomes: list[str], filename: str | Path) -> None:
    """Write to an OrthoFinder-style orthogroups.tsv file.

    This function takes an iterator of orthogroup results and saves them to a tab-separated values (TSV)
    file. The first column contains the SOG ID, followed by columns for each genome in the `all_genomes` list.
    For each genome, the function lists the gene IDs associated with the SOG.

    Args:
        results_gen: An iterator that yields tuples containing an orthogroup ID and a dictionary mapping
                     genomes to lists of genes.
        all_genomes (list[str]): A list of genome names corresponding to the columns in the TSV file.
        filename (str | Path): The path where the resulting TSV file will be saved.

    Returns:
        None: This function does not return any value; it writes directly to the specified file.
    """
    with open(filename, "w", encoding="utf-8") as f:
        # Write Header
        header = ["SOG_ID"] + all_genomes
        f.write("\t".join(header) + "\n")

        # Write Rows
        for sog in results_gen:
            row = [sog.id]
            sog_dict = defaultdict(list)
            for g in sog:
                sog_dict[g.genome.name].append(g)
            for g_name in all_genomes:
                matching_genes = sog_dict.get(g_name, None)
                row.append(",".join(g.id for g in matching_genes) if matching_genes else "")
            f.write("\t".join(row) + "\n")
