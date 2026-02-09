#!/usr/bin/env python3
"""
orthoSynAssign - Ortholog Synteny Assignment Tool

Main program for analyzing orthologous groups and synteny using OrthoFinder
results and genome annotation files (GFF3/GTF).
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import TYPE_CHECKING

from .lib.parsers import BedParser, GFFParser, GTFParser, read_orthofinder_table, save_results_tsv


if TYPE_CHECKING:
    from .lib.parsers import AnnotationParser


def setup_logging(verbose: bool = False):
    """Configure logging for the application."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Ortholog Synteny Assignment Tool - Analyze orthologous groups using "
        "OrthoFinder results and genome annotations",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Using GFF3 files
  orthoSynAssign.py --gff annotations/ --orthofinder Orthogroups.tsv

  # Using GTF files
  orthoSynAssign.py --gtf gtf_files/ --orthofinder Orthogroups.tsv

  # With verbose output
  orthoSynAssign.py --gff annotations/ --orthofinder Orthogroups.tsv -v
        """,
    )

    # OrthoFinder input
    parser.add_argument(
        "-i",
        "--og_input",
        dest="og_input",
        type=Path,
        required=True,
        help="Path to OrthoFinder Orthogroups.tsv file",
    )

    # Input format group (mutually exclusive)
    format_group = parser.add_mutually_exclusive_group(required=True)
    format_group.add_argument(
        "--bed",
        type=Path,
        metavar=("file", "files"),
        nargs="+",
        help="Path of BED formatted genome annotation files",
    )
    format_group.add_argument(
        "--gff",
        type=Path,
        metavar=("file", "files"),
        nargs="+",
        help="Path of GFF3 formatted genome annotation files",
    )
    format_group.add_argument(
        "--gtf",
        type=Path,
        metavar=("file", "files"),
        nargs="+",
        help="Path of GTF formatted genome annotation files",
    )

    parser.add_argument(
        "-w",
        "--window",
        dest="window",
        type=int,
        default=8,
        help="Controls how many total genes are considered when determining synteny for a single gene",
    )

    parser.add_argument(
        "-r",
        "--ratio_threshold",
        dest="ratio_threshold",
        type=float,
        default=0.5,
        help="Controls how many genes within a window must provide synteny support to classify the genes being compared as syntenous",
    )

    # Optional arguments
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        type=Path,
        default="OrthoRefined_SOGs.tsv",
        help="Output of results",
    )

    parser.add_argument(
        "--skip_single_orthologs",
        dest="skip_single_orthologs",
        action="store_true",
        help="Skip single orthologs",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")

    return parser.parse_args()


def _validate_annotations(args, logger) -> list[AnnotationParser]:
    """
    Validate input files and directories.

    Args:
        args: Parsed command line arguments
    """
    # Check annotation files
    if args.bed:
        files = args.bed
        parser = BedParser
    elif args.gff:
        files = args.gff
        parser = GFFParser
    else:
        files = args.gtf
        parser = GTFParser

    annotations = []
    for file in files:
        logger.debug(f"Reading annotation from: {file}")
        annotations.append(parser(file))

    return annotations


def _validate_orthogroup(args) -> Path:
    # Check OrthoFinder file
    file = Path(args.og_input)
    if not file.exists():
        raise FileNotFoundError(f"Orthogroup file not found: {file}")
    if not file.is_file():
        raise ValueError(f"Orthogroup file must be a regular file: {file}")
    return file


def main():
    """Main entry point for orthoSynAssign."""
    # Parse arguments
    args = parse_arguments()

    # Setup logging
    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)

    logger.info("Starting orthoSynAssign")
    logger.info(f"Command: {' '.join(sys.argv)}")

    try:
        # Validate inputs
        annotations = _validate_annotations(args, logger)
        og_file = _validate_orthogroup(args)

        # Read gff
        genomes = {}
        for annotation in annotations:
            genome = annotation.parse()
            genomes[genome.sample_name] = genome

        # Read OrthoFinder orthogroups
        logger.info(f"Reading OrthoFinder data from: {og_file}")
        orthogroups = read_orthofinder_table(og_file, genomes)
        logger.info(f"Loaded {len(orthogroups)} orthogroups")

        # # Create output directory
        # output_dir = Path(args.output)
        # output_dir.mkdir(parents=True, exist_ok=True)
        # logger.info(f"Output directory: {output_dir}")

        # Perform synteny analysis
        logger.info("Analysis placeholder - implement your analysis logic here")

        final_sog_results = []
        global_sog_counter = 1

        for og in orthogroups:
            # Each 'og' is an instance of the Orthogroup class we built
            # It handles its own refinement logic internally
            refined_sogs = og.get_refined_sogs(
                window_size=args.window, ratio_threshold=args.ratio_threshold, skip_single_ortholog=args.skip_single_orthologs
            )

            # If refinement results in empty (no syntenic support), skip it
            if not refined_sogs:
                continue

            # Store results with a sub-identifier (HOG0001.0, HOG0001.1, etc.)
            for sog_dict in refined_sogs:
                # Format as SOG000001, SOG000002, etc.
                sog_id = f"SOG{global_sog_counter:06d}"

                # Optional: Keep the original HOG ID as metadata in the tuple
                final_sog_results.append((f"{sog_id}.{og.og_id}", sog_dict))
                global_sog_counter += 1

        # --- 4. OUTPUT TO TSV ---
        save_results_tsv(final_sog_results, args.output)
        logger.info(f"Refinement complete. Results saved to {args.output}")

        logger.info("orthoSynAssign completed successfully")

    except Exception as e:
        logger.exception(f"An error occurred: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
