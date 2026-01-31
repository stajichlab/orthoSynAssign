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

from .lib import read_orthofinder_table
from .lib.parsers import read_gff_folder, read_gtf_folder, save_results_tsv


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

    # OrthoFinder input
    parser.add_argument(
        "-i--og_input",
        dest="og_input",
        type=str,
        required=True,
        help="Path to OrthoFinder Orthogroups.tsv file",
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
        type=str,
        default="output",
        help="Output directory for results (default: output)",
    )

    parser.add_argument(
        "--skip_single_orthologs",
        dest="skip_single_orthologs",
        action="store_true",
        help="Output directory for results (default: output)",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")

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
        errors.append(f"OrthoFinder file not found: {args.og_input}")

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
        logger.info(f"Reading OrthoFinder data from: {args.og_input}")
        orthogroups = read_orthofinder_table(args.og_input)
        logger.info(f"Loaded {len(orthogroups)} orthogroups")

        # Create output directory
        output_dir = Path(args.output)
        output_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Output directory: {output_dir}")

        # TODO: Implement analysis logic here
        # This is where you would add your synteny analysis, comparison, etc.
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
                final_sog_results.append((sog_id, og.og_id, sog_dict))
                global_sog_counter += 1

        # --- 4. OUTPUT TO TSV ---
        output_file = "OrthoRefined_SOGs.tsv"
        save_results_tsv(final_sog_results, output_file)
        logger.info(f"Refinement complete. Results saved to {output_file}")

        logger.info("orthoSynAssign completed successfully")

    except Exception as e:
        logger.exception(f"An error occurred: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
