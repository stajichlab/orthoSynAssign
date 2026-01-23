"""
File parsers for GTF, GFF3, and OrthoFinder ortholog tables.
"""

import os
from pathlib import Path
from typing import List, Dict, Any, Optional
import logging
import gzip

logger = logging.getLogger(__name__)


class GFFFeature:
    """Represents a feature from a GFF3/GTF file."""

    def __init__(
        self,
        seqid: str,
        source: str,
        feature_type: str,
        start: int,
        end: int,
        score: Optional[str],
        strand: str,
        phase: str,
        attributes: Dict[str, str],
    ):
        self.seqid = seqid
        self.source = source
        self.type = feature_type
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.phase = phase
        self.attributes = attributes

    def __repr__(self):
        return (
            f"GFFFeature({self.seqid}:{self.start}-{self.end} "
            f"{self.strand} {self.type})"
        )


def parse_gff3_attributes(attr_string: str) -> Dict[str, str]:
    """
    Parse GFF3 attributes field.

    GFF3 format: ID=gene1;Name=MyGene;Note=Some note

    Args:
        attr_string: The attributes column from a GFF3 line

    Returns:
        Dictionary of attribute key-value pairs
    """
    attributes = {}
    if not attr_string or attr_string == ".":
        return attributes

    for item in attr_string.split(";"):
        item = item.strip()
        if "=" in item:
            key, value = item.split("=", 1)
            attributes[key] = value

    return attributes


def parse_gtf_attributes(attr_string: str) -> Dict[str, str]:
    """
    Parse GTF attributes field.

    GTF format: gene_id "ENSG00000123"; transcript_id "ENST00000456";

    Args:
        attr_string: The attributes column from a GTF line

    Returns:
        Dictionary of attribute key-value pairs
    """
    attributes = {}
    if not attr_string or attr_string == ".":
        return attributes

    # Split by semicolon and process each key-value pair
    for item in attr_string.split(";"):
        item = item.strip()
        if item:
            parts = item.split(None, 1)  # Split on first whitespace
            if len(parts) == 2:
                key = parts[0]
                value = parts[1].strip('"')  # Remove quotes
                attributes[key] = value

    return attributes


def read_gff3(file_path: str) -> List[GFFFeature]:
    """
    Read a GFF3 formatted genome annotation file.

    Args:
        file_path: Path to the GFF3 file

    Returns:
        List of GFFFeature objects

    Raises:
        FileNotFoundError: If the file doesn't exist
        ValueError: If the file format is invalid
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"GFF3 file not found: {file_path}")

    features = []
    logger.info(f"Reading GFF3 file: {file_path}")

    if file_path.endswith(".gz"):
        f = gzip.open(file_path, "rt", encoding="utf-8")
    else:
        f = open(file_path, "r", encoding="utf-8")
    with f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()

            # Skip empty lines and comments
            if not line or line.startswith("#"):
                continue

            # Parse the line
            fields = line.split("\t")
            if len(fields) != 9:
                logger.warning(f"Line {line_num}: Expected 9 fields, got {len(fields)}")
                continue

            (
                seqid,
                source,
                feature_type,
                start,
                end,
                score,
                strand,
                phase,
                attributes,
            ) = fields

            try:
                feature = GFFFeature(
                    seqid=seqid,
                    source=source,
                    feature_type=feature_type,
                    start=int(start),
                    end=int(end),
                    score=score if score != "." else None,
                    strand=strand,
                    phase=phase,
                    attributes=parse_gff3_attributes(attributes),
                )
                features.append(feature)
            except ValueError as e:
                logger.warning(f"Line {line_num}: Error parsing coordinates - {e}")
                continue

    logger.info(f"Loaded {len(features)} features from {file_path}")
    return features


def read_gtf(file_path: str) -> List[GFFFeature]:
    """
    Read a GTF formatted genome annotation file.

    Args:
        file_path: Path to the GTF file

    Returns:
        List of GFFFeature objects

    Raises:
        FileNotFoundError: If the file doesn't exist
        ValueError: If the file format is invalid
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"GTF file not found: {file_path}")

    features = []
    logger.info(f"Reading GTF file: {file_path}")

    if file_path.endswith(".gz"):
        f = gzip.open(file_path, "rt", encoding="utf-8")
    else:
        f = open(file_path, "r", encoding="utf-8")
    with f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()

            # Skip empty lines and comments
            if not line or line.startswith("#"):
                continue

            # Parse the line
            fields = line.split("\t")
            if len(fields) != 9:
                logger.warning(f"Line {line_num}: Expected 9 fields, got {len(fields)}")
                continue

            (
                seqid,
                source,
                feature_type,
                start,
                end,
                score,
                strand,
                phase,
                attributes,
            ) = fields

            try:
                feature = GFFFeature(
                    seqid=seqid,
                    source=source,
                    feature_type=feature_type,
                    start=int(start),
                    end=int(end),
                    score=score if score != "." else None,
                    strand=strand,
                    phase=phase,
                    attributes=parse_gtf_attributes(attributes),
                )
                features.append(feature)
            except ValueError as e:
                logger.warning(f"Line {line_num}: Error parsing coordinates - {e}")
                continue

    logger.info(f"Loaded {len(features)} features from {file_path}")
    return features


def read_orthofinder_table(file_path: str) -> Dict[str, Dict[str, List[str]]]:
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
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"OrthoFinder file not found: {file_path}")

    orthogroups = {}
    logger.info(f"Reading OrthoFinder table: {file_path}")

    if file_path.endswith(".gz"):
        f = gzip.open(file_path, "rt", encoding="utf-8")
    else:
        f = open(file_path, "r", encoding="utf-8")
    with f:
        # Read header to get species names
        header = f.readline().strip().split("\t")
        if len(header) < 2:
            raise ValueError("Invalid OrthoFinder file: expected at least 2 columns")

        # First column is "Orthogroup", rest are species names
        species_names = header[1:]
        logger.info(f"Found {len(species_names)} species: {', '.join(species_names)}")

        # Read orthogroup data
        for line_num, line in enumerate(f, 2):  # Start at 2 since we read header
            line = line.strip("\n")
            if not line:
                continue

            fields = line.split("\t")
            if len(fields) != len(header):
                logger.warning(
                    f"Line {line_num}: Expected {len(header)} fields, got {len(fields)}"
                )
                continue

            orthogroup_id = fields[0]
            orthogroups[orthogroup_id] = {}

            # Parse genes for each species
            for i, species in enumerate(species_names, 1):
                genes_str = fields[i].strip()
                if genes_str:
                    # Genes are usually separated by commas or spaces
                    genes = [
                        g.strip()
                        for g in genes_str.replace(",", " ").split()
                        if g.strip()
                    ]
                    orthogroups[orthogroup_id][species] = genes
                else:
                    orthogroups[orthogroup_id][species] = []

    logger.info(f"Loaded {len(orthogroups)} orthogroups from {file_path}")
    return orthogroups


def read_gff_folder(
    folder_path: str, file_extension: str = ".gff"
) -> Dict[str, List[GFFFeature]]:
    """
    Read all GFF3 files from a folder.

    Args:
        folder_path: Path to folder containing GFF3 files
        file_extension: File extension to filter (default: .gff3)

    Returns:
        Dictionary mapping file names (without extension) to lists of features
    """
    if not os.path.isdir(folder_path):
        raise NotADirectoryError(f"Not a directory: {folder_path}")

    results = {}
    folder = Path(folder_path)
    file_paths = []
    for file_path in folder.glob(f"*{file_extension}*"):
        if file_path.is_file() and (
            file_path.name.endswith(file_extension)
            or file_path.name.endswith(f"{file_extension}3")
            or file_path.name.endswith(f"{file_extension}3.gz")
            or file_path.name.endswith(f"{file_extension}.gz")
        ):
            species_name = file_path.stem
            results[species_name] = read_gff3(str(file_path))

    logger.info(f"Loaded annotations for {len(results)} species from {folder_path}")
    return results


def read_gtf_folder(
    folder_path: str, file_extension: str = ".gtf"
) -> Dict[str, List[GFFFeature]]:
    """
    Read all GTF files from a folder.

    Args:
        folder_path: Path to folder containing GTF files
        file_extension: File extension to filter (default: .gtf)

    Returns:
        Dictionary mapping file names (without extension) to lists of features
    """
    if not os.path.isdir(folder_path):
        raise NotADirectoryError(f"Not a directory: {folder_path}")

    results = {}
    folder = Path(folder_path)

    for file_path in folder.glob(f"*{file_extension}*"):
        if file_path.is_file() and (
            file_path.name.endswith(file_extension)
            or file_path.name.endswith(f"{file_extension}.gz")
        ):
            species_name = file_path.stem
            logger.info(f"Reading {species_name} from {file_path}")
            results[species_name] = read_gtf(str(file_path))

    logger.info(f"Loaded annotations for {len(results)} species from {folder_path}")
    return results
