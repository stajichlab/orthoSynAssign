#!/usr/bin/env python3
"""
orthoSynAssign - Ortholog Synteny Assignment Tool

Main program for analyzing orthologous groups and synteny using OrthoFinder
results and genome annotation files (GFF3/GTF).
"""

import argparse
import logging
import sys
from pathlib import Path

from orthoSynAssign.lib import read_gtf, read_gff3, read_orthofinder_table
from orthoSynAssign.lib.parsers import read_gff_folder, read_gtf_folder


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
  orthoSynAssign.py --gff_folder annotations/ --orthofinder Orthogroups.tsv

  # Using GTF files
  orthoSynAssign.py --gtf_folder gtf_files/ --orthofinder Orthogroups.tsv

  # With verbose output
  orthoSynAssign.py --gff_folder annotations/ --orthofinder Orthogroups.tsv -v
        """,
    )

    # Input format group (mutually exclusive)
    format_group = parser.add_mutually_exclusive_group(required=True)
    format_group.add_argument(
        "--gff_folder",
        type=str,
        help="Path to folder containing GFF3 formatted genome annotation files",
    )
    format_group.add_argument(
        "--gtf_folder",
        type=str,
        help="Path to folder containing GTF formatted genome annotation files",
    )

    # OrthoFinder input
    parser.add_argument(
        "--orthofinder",
        type=str,
        required=True,
        help="Path to OrthoFinder Orthogroups.tsv file",
    )

    # Optional arguments
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="output",
        help="Output directory for results (default: output)",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose logging"
    )

    return parser.parse_args()


def validate_inputs(args):
    """
    Validate input files and directories.

    Args:
        args: Parsed command line arguments

    Returns:
        True if all inputs are valid, False otherwise
    """
    errors = []

    # Check annotation folder
    if args.gff_folder:
        if not Path(args.gff_folder).is_dir():
            errors.append(f"GFF folder not found: {args.gff_folder}")
    elif args.gtf_folder:
        if not Path(args.gtf_folder).is_dir():
            errors.append(f"GTF folder not found: {args.gtf_folder}")

    # Check OrthoFinder file
    if not Path(args.orthofinder).is_file():
        errors.append(f"OrthoFinder file not found: {args.orthofinder}")

    if errors:
        for error in errors:
            logging.error(error)
        return False

    return True


def main():
    """Main entry point for orthoSynAssign."""
    # Parse arguments
    args = parse_arguments()

    # Setup logging
    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)

    logger.info("Starting orthoSynAssign")
    logger.info(f"Command: {' '.join(sys.argv)}")

    # Validate inputs
    if not validate_inputs(args):
        logger.error("Input validation failed. Exiting.")
        sys.exit(1)

    try:
        # Read annotation files
        annotations = {}
        if args.gff_folder:
            logger.info(f"Reading GFF3 files from: {args.gff_folder}")
            annotations = read_gff_folder(args.gff_folder)
        elif args.gtf_folder:
            logger.info(f"Reading GTF files from: {args.gtf_folder}")
            annotations = read_gtf_folder(args.gtf_folder)

        logger.info(f"Loaded annotations for {len(annotations)} species")
        for species, features in annotations.items():
            logger.info(f"  {species}: {len(features)} features")

        # Read OrthoFinder orthogroups
        logger.info(f"Reading OrthoFinder data from: {args.orthofinder}")
        orthogroups = read_orthofinder_table(args.orthofinder)
        logger.info(f"Loaded {len(orthogroups)} orthogroups")

        # Create output directory
        output_dir = Path(args.output)
        output_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Output directory: {output_dir}")

        # TODO: Implement analysis logic here
        # This is where you would add your synteny analysis, comparison, etc.
        logger.info("Analysis placeholder - implement your analysis logic here")

        # Example: Print summary statistics
        total_genes_in_orthogroups = sum(
            len(genes) for og in orthogroups.values() for genes in og.values()
        )
        logger.info(f"Total genes in orthogroups: {total_genes_in_orthogroups}")

        logger.info("orthoSynAssign completed successfully")

    except Exception as e:
        logger.exception(f"An error occurred: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
