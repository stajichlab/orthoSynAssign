#!/usr/bin/env python3
"""Convert GFF3 or GTF to BED format, mapping proteins to their genes."""

import argparse
import gzip
import re
import sys
from collections import defaultdict


def get_attr(attributes_str, key):
    """Extract a value from a GFF3 attributes field (key=value;...)."""
    for attr in attributes_str.split(";"):
        if attr.startswith(key + "="):
            return attr.split("=", 1)[1]
    return ""


def get_gtf_attr(attributes_str, key):
    """Extract a value from a GTF attributes field (key "value"; ...)."""
    for attr in attributes_str.strip().split(";"):
        attr = attr.strip()
        if attr.startswith(key + " "):
            return attr[len(key) + 1 :].strip().strip('"')
    return ""


# Predefined format styles.
#
# Each style declares the same fields so a single algorithm per dialect can
# handle all of them:
#   dialect          "gff3" (ID=/Parent= attributes) or "gtf" (key "value")
#   gene_type        feature type that defines a gene (and its coordinates)
#   protein_type     feature type that carries protein identity
#   protein_id_attr  attribute holding the protein id; for gff3, None means
#                    use the protein feature's own ID
#   gene_id_attr     (gtf only) attribute holding the gene id
FORMATS = {
    "eupathdb": {
        "description": "EuPathDB GFF3 (e.g., PlasmoDB, TriTrypDB)",
        "dialect": "gff3",
        "gene_type": "protein_coding_gene",
        "protein_type": "CDS",
        "protein_id_attr": "protein_source_id",
    },
    "sequence_ontology": {
        "description": "GFF3-SequenceOntology (e.g., funannotate)",
        "dialect": "gff3",
        "gene_type": "gene",
        "protein_type": "mRNA",
        "protein_id_attr": None,
    },
    "refseq": {
        "description": "NCBI RefSeq GFF3 (protein accessions from CDS protein_id)",
        "dialect": "gff3",
        "gene_type": "gene",
        "protein_type": "CDS",
        "protein_id_attr": "protein_id",
    },
    "ensembl": {
        "description": "Ensembl/EnsemblGenomes GFF3 (protein ids from CDS protein_id)",
        "dialect": "gff3",
        "gene_type": "gene",
        "protein_type": "CDS",
        "protein_id_attr": "protein_id",
    },
    "gtf": {
        "description": "GTF (RefSeq/Ensembl; gene_id and protein_id on each row)",
        "dialect": "gtf",
        "gene_type": "gene",
        "protein_type": "CDS",
        "gene_id_attr": "gene_id",
        "protein_id_attr": "protein_id",
    },
}

# Default style to use for each detected dialect when --style is not given.
DEFAULT_STYLE = {
    "gff3": "sequence_ontology",
    "gtf": "gtf",
}


def parse_args():
    parser = argparse.ArgumentParser(description="Convert GFF3 or GTF to BED format, mapping proteins to their genes.")
    parser.add_argument("gff_file", nargs="?", help="GFF3 or GTF file (gzipped or plain text)")
    parser.add_argument(
        "--style", choices=FORMATS.keys(), default=None, help="Format style (default: auto-selected from the detected dialect)"
    )
    parser.add_argument(
        "--format", choices=["auto", "gff3", "gtf"], default="auto", help="Input dialect (default: auto-detect from file content)"
    )
    parser.add_argument("--list-styles", action="store_true", help="List available format styles and exit")
    parser.add_argument(
        "--use-gene-id", action="store_true", help="Use gene ID as name column instead of transcript ID (funannotate only)"
    )
    parser.add_argument(
        "--collapse-transcripts",
        action="store_true",
        help="Collapse all transcripts per gene into one line with semicolon-separated IDs (funannotate only)",
    )

    return parser.parse_args()


def list_styles():
    print("Available format styles:\n")
    for name, config in FORMATS.items():
        print(f"  {name:20} [{config['dialect']:4}] {config['description']}")
    print()


def open_gff(gff_file):
    """Open a GFF3/GTF file, transparently handling gzip compression."""
    with open(gff_file, "rb") as probe:
        is_gzip = probe.read(2) == b"\x1f\x8b"

    if is_gzip:
        return gzip.open(gff_file, "rt")
    return open(gff_file, "r")


def detect_dialect(gff_file):
    """Sniff the attribute syntax of the first data line to pick a dialect."""
    with open_gff(gff_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            attributes = fields[8]
            # GTF uses `key "value"`; GFF3 uses `key=value`.
            if re.search(r'\w+\s+"', attributes):
                return "gtf"
            if "=" in attributes:
                return "gff3"
    return "gff3"


def gff_to_gene_proteins(lines, config):
    """Map genes to their protein ids from GFF3 lines.

    Returns a tuple of:
      gene_coords      gene_id -> (chrom, start, end) with 0-based BED start
      gene_to_proteins gene_id -> list of protein ids (ordered by appearance)
    """
    gene_type = config["gene_type"]
    protein_type = config["protein_type"]
    protein_id_attr = config["protein_id_attr"]

    gene_coords = {}  # gene_id -> (chrom, start, end)
    parent_map = {}  # feature_id -> parent_id
    protein_features = []  # [(protein_id, parent_id, order)]
    gene_to_proteins = defaultdict(list)  # gene_id -> list of protein ids
    feature_order = 0

    # PASS 1: read features, record gene coordinates, parent links, and the
    # features that carry protein identity.
    for line in lines:
        if line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue

        chrom, source, ftype, start, end, score, strand, phase, attributes = fields

        feature_id = get_attr(attributes, "ID")
        parent = get_attr(attributes, "Parent")

        if feature_id:
            # Store gene coordinates (converting to 0-based BED format).
            if ftype == gene_type:
                gene_coords[feature_id] = (chrom, int(start) - 1, int(end))

            if parent:
                parent_map[feature_id] = parent

        # Record protein-bearing features. The protein id comes from a named
        # attribute, or from the feature's own ID when no attribute is set.
        if ftype == protein_type:
            protein_id = get_attr(attributes, protein_id_attr) if protein_id_attr else feature_id
            if protein_id and parent:
                protein_features.append((protein_id, parent, feature_order))
                feature_order += 1

    # PASS 2: climb the parent chain from each protein feature up to its gene.
    for protein_id, parent, order in protein_features:
        current = parent
        while current and current not in gene_coords:
            current = parent_map.get(current, "")

        if current in gene_coords:
            gene_to_proteins[current].append((protein_id, order))

    # Sort by appearance order and keep only IDs
    for gene_id in gene_to_proteins:
        gene_to_proteins[gene_id].sort(key=lambda x: x[1])
        gene_to_proteins[gene_id] = list(dict.fromkeys([pid for pid, _ in gene_to_proteins[gene_id]]))

    return gene_coords, gene_to_proteins


def gtf_to_gene_proteins(lines, config):
    """Map genes to their protein ids from GTF lines.

    GTF carries gene_id (and protein_id) directly on every row, so no parent
    climbing is needed. Returns the same tuple as gff_to_gene_proteins.
    """
    gene_type = config["gene_type"]
    protein_type = config["protein_type"]
    gene_id_attr = config["gene_id_attr"]
    protein_id_attr = config["protein_id_attr"]

    gene_coords = {}  # gene_id -> (chrom, start, end)
    gene_to_proteins = defaultdict(set)  # gene_id -> set of protein ids

    for line in lines:
        if line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue

        chrom, source, ftype, start, end, score, strand, phase, attributes = fields

        if ftype == gene_type:
            gene_id = get_gtf_attr(attributes, gene_id_attr)
            if gene_id:
                gene_coords[gene_id] = (chrom, int(start) - 1, int(end))
        elif ftype == protein_type:
            gene_id = get_gtf_attr(attributes, gene_id_attr)
            protein_id = get_gtf_attr(attributes, protein_id_attr)
            if gene_id and protein_id:
                gene_to_proteins[gene_id].add(protein_id)

    # Keep only proteins whose gene has coordinates (matches GFF3 behavior).
    resolved = defaultdict(set)
    for gene_id, proteins in gene_to_proteins.items():
        if gene_id in gene_coords:
            resolved[gene_id] = proteins

    return gene_coords, resolved


def main():
    args = parse_args()

    if args.list_styles:
        list_styles()
        sys.exit(0)

    if not args.gff_file:
        sys.exit("error: gff_file is required (or use --list-styles)")

    dialect = args.format if args.format != "auto" else detect_dialect(args.gff_file)

    if args.style:
        config = FORMATS[args.style]
        if config["dialect"] != dialect:
            sys.exit(
                f"error: style '{args.style}' is {config['dialect']} but the "
                f"input looks like {dialect}; pass --format to override"
            )
    else:
        config = FORMATS[DEFAULT_STYLE[dialect]]

    parse_dialect = gtf_to_gene_proteins if dialect == "gtf" else gff_to_gene_proteins

    with open_gff(args.gff_file) as f:
        gene_coords, gene_to_proteins = parse_dialect(f, config)

    # Build BED records.
    bed_records = []
    for gene_id, proteins in gene_to_proteins.items():
        chrom, start, end = gene_coords[gene_id]

        # For funannotate (sequence_ontology), use transcript ID by default
        if config.get("protein_type") == "mRNA" and (args.style == "sequence_ontology" or (not args.style and dialect == "gff3")):
            if args.use_gene_id:
                # Use gene ID as the name
                name = gene_id
            elif args.collapse_transcripts:
                # Use semicolon-separated transcript IDs
                name = ";".join(proteins)
            else:
                # Use only first transcript ID
                name = proteins[0] if proteins else gene_id
            bed_records.append((chrom, start, end, name))
        else:
            # Legacy behavior for other formats: use gene_id and list all proteins
            protein_list = ";".join(sorted(proteins))
            bed_records.append((chrom, start, end, gene_id, protein_list))

    # Sort by chromosome and position; end and name break ties so the
    # output is deterministic.
    if len(bed_records) > 0 and len(bed_records[0]) == 4:
        # 4-column BED format
        bed_records.sort(key=lambda x: (x[0], x[1], x[2], x[3]))
        for chrom, start, end, name in bed_records:
            print(f"{chrom}\t{start}\t{end}\t{name}")
    else:
        # 5-column format (legacy)
        bed_records.sort(key=lambda x: (x[0], x[1], x[2], x[3]))
        for chrom, start, end, gene_id, proteins in bed_records:
            print(f"{chrom}\t{start}\t{end}\t{gene_id}\t{proteins}")


if __name__ == "__main__":
    main()
