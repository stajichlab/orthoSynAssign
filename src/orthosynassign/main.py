#!/usr/bin/env python3
"""
The orthosynassign CLI entry point.
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
import textwrap
import time
from pathlib import Path
from typing import TYPE_CHECKING

from . import AUTHOR, VERSION
from . import __doc__ as _module_doc
from .lib.parsers import BedParser, read_orthofinder_table, save_results_tsv

if TYPE_CHECKING:
    from .lib.parsers import AnnotationParser

_EPILOG = textwrap.dedent(f"""\
Examples:

# Specify bed files separately:
orthosynassign -i orthogroup.tsv --bed file1.bed file2.bed file3.bed

# Specify all bed files in a directory:
orthosynassign -i orthogroup.tsv --bed *.bed

# Specify output file name for results:
orthosynassign -i orthogroup.tsv --bed *.bed -o Refined_orthologs.tsv

# Specify window size and ratio threshold:
orthosynassign -i orthogroup.tsv --bed *.bed -w 10 -r 0.8

# With verbose output:
orthosynassign -i orthogroup.tsv --bed *.bed -v

Written by {AUTHOR}
""")


def main():
    """Main entry point for orthoSynAssign."""
    # Parse arguments
    args = _parse_arguments()

    # Setup logging
    _setup_logging(args.verbose)
    logger = logging.getLogger(__name__)

    logger.info("Starting orthoSynAssign")
    logger.debug(f"Command: {' '.join(sys.argv)}")

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
        logger.debug(f"Loaded {len(orthogroups)} orthogroups")

        # Create output directory
        output_dir = Path(args.output).parent
        output_dir.mkdir(parents=True, exist_ok=True)
        logger.debug(f"Output directory: {output_dir}")

        # Perform synteny analysis
        logger.info("Refining orthogroups by pairwise synteny analysis.")

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

            # Store results with an identifier
            for sog_dict in refined_sogs:
                # Format as SOG000001, SOG000002, etc.
                sog_id = f"SOG{global_sog_counter:06d}"

                final_sog_results.append((f"{sog_id}.{og.og_id}", sog_dict))
                global_sog_counter += 1

        save_results_tsv(final_sog_results, args.output)
        logger.info(f"Refinement complete. Results saved to {args.output}")

        logger.info("orthoSynAssign completed successfully")

    except Exception as e:
        logger.error(f"An error occurred: {e}")
        sys.exit(1)


class _CustomHelpFormatter(argparse.RawDescriptionHelpFormatter):
    """
    Custom help formatter for argparse to enhance text wrapping and default value display.

    This formatter adjusts the text wrapping for better readability and ensures that the default value of an argument is displayed
    in the help message if not already present.
    """

    def _get_help_string(self, action):
        """Allow additional message after default parameter displayed."""
        help = action.help
        pattern = r"\(default: .+\)"
        if re.search(pattern, action.help) is None:
            if action.default not in [argparse.SUPPRESS, None, False]:
                defaulting_nargs = [argparse.OPTIONAL, argparse.ZERO_OR_MORE]
                if action.option_strings or action.nargs in defaulting_nargs:
                    help += " (default: %(default)s)"
        return help


def _parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description=_module_doc,
        formatter_class=_CustomHelpFormatter,
        epilog=_EPILOG,
        add_help=False,
    )
    req_args = parser.add_argument_group("Required arguments")
    # OrthoFinder input
    req_args.add_argument(
        "-i",
        "--og_input",
        dest="og_input",
        type=Path,
        required=True,
        help="Path to OrthoFinder Orthogroups.tsv file",
    )

    # Input format group (mutually exclusive)
    req_args.add_argument(
        "--bed",
        type=Path,
        required=True,
        metavar=("file", "files"),
        nargs="+",
        help="Path of BED formatted genome annotation files",
    )

    opt_args = parser.add_argument_group("Options")
    opt_args.add_argument(
        "-w",
        "--window",
        dest="window",
        type=int,
        default=8,
        help="Controls how many total genes are considered when determining synteny for a single gene",
    )

    opt_args.add_argument(
        "-r",
        "--ratio_threshold",
        dest="ratio_threshold",
        type=float,
        default=0.5,
        help=textwrap.dedent("""
            Controls how many genes within a window must provide synteny support
            to classify the genes being compared as syntenous
        """),
    )

    opt_args.add_argument(
        "-o",
        "--output",
        dest="output",
        type=Path,
        default="Refined_SOGs-%s.tsv" % time.strftime("%Y%m%d-%H%M%S", time.gmtime()),
        help="Output of results (default: Refined_SOGs-[YYYYMMDD-HHMMSS].tsv (UTC timestamp))",
    )

    # opt_args.add_argument(
    #     "--skip_single_orthologs",
    #     dest="skip_single_orthologs",
    #     action="store_true",
    #     help="Skip single orthologs",
    # )
    opt_args.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    opt_args.add_argument("-V", "--version", action="version", version=VERSION)
    opt_args.add_argument("-h", "--help", action="help", help="show this help message and exit")

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


def _setup_logging(verbose: bool = False):
    """Configure logging for the application."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


if __name__ == "__main__":
    main()
