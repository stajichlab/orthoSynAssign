#!/usr/bin/env python3
"""
The orthosynassign CLI entry point.
"""

from __future__ import annotations

import argparse
import logging
import multiprocessing
import os
import re
import sys
import textwrap
import time
import traceback
from functools import partial
from pathlib import Path
from typing import TYPE_CHECKING, Iterator

from tqdm import tqdm

from . import AUTHOR, VERSION
from . import __doc__ as _module_doc
from .lib.parsers import BedParser, read_orthofinder_table, save_results_tsv

if TYPE_CHECKING:
    from .lib._gene import Gene, Genome
    from .lib._orthogroup import Orthogroup
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

_WORKER_OGS: list[Orthogroup] = []
_WORKER_GENOME_MAP: dict[str, Genome] = {}
_AVAIL_CPUS = int(os.environ.get("SLURM_CPUS_ON_NODE", os.cpu_count()))


def main(args: list[str] | None = None) -> int:
    """Main entry point for orthoSynAssign."""
    # Parse arguments
    args = _parse_arguments() or sys.argv[1:]

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
        results_stream = _generate_sog_results(orthogroups, genomes, args, cpus=args.threads)

        save_results_tsv(results_stream, list(genomes.keys()), args.output)
        logger.info(f"Refinement complete. Results saved to {args.output}")

        logger.info("orthoSynAssign completed successfully")

    except KeyboardInterrupt:
        logger.warning("Terminated by user.")
        return 1

    except Exception as e:
        logger.error(f"An error occurred: {e}")
        logger.debug("%s", traceback.format_exc())
        return 1

    return 0


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

    opt_args.add_argument(
        "--skip_single_orthologs",
        dest="skip_single_orthologs",
        action="store_true",
        help="Skip single orthologs",
    )
    opt_args.add_argument("-t", "--threads", type=int, default=min(_AVAIL_CPUS, 4), help="Number of cpus to use")
    opt_args.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    opt_args.add_argument("-V", "--version", action="version", version=VERSION)
    opt_args.add_argument("-h", "--help", action="help", help="show this help message and exit")

    return parser.parse_args()


def _setup_logging(verbose: bool = False):
    """Configure logging for the application."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


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


def _generate_sog_results(
    orthogroups: list[Orthogroup], genome_data: dict[str, Genome], args, *, cpus: int = 1
) -> Iterator[tuple[str, dict[Genome, list[Gene]]]]:
    """
    Logic-only: Processes orthogroups and yields results one by one.
    """
    global_sog_counter = 1
    total_ogs = len(orthogroups)
    indices = list(range(total_ogs))
    desc = f"Processing with {cpus} cpu{'s' if cpus > 1 else ''}"

    if cpus == 1:
        _init_worker(orthogroups, genome_data)

        for index in tqdm(indices, desc=desc, unit="og"):
            _, refined_sogs = _process_og_task(index, args)

            for sog_dict in refined_sogs:
                sog_id = f"SOG{global_sog_counter:06d}.{orthogroups[index].og_id}"
                yield sog_id, sog_dict
                global_sog_counter += 1
        return

    # Use 'spawn' to ensure Windows compatibility
    ctx = multiprocessing.get_context("spawn")

    opt_chunksize = _calculate_optimal_chunksize(total_ogs, cpus)
    # Start the Pool with the shared memory pointer
    pool = ctx.Pool(processes=cpus, initializer=_init_worker, initargs=(orthogroups, genome_data))

    try:
        # Pass the (og, args) tuple to imap_unordered
        worker_func = partial(_process_og_task, args=args)

        pbar = tqdm(pool.imap(worker_func, indices, chunksize=opt_chunksize), total=total_ogs, desc=desc, unit="og")

        for og_id, refined_sogs in pbar:
            for sog_dict in refined_sogs:
                sog_id_str = f"SOG{global_sog_counter:06d}.{og_id}"
                yield sog_id_str, sog_dict
                global_sog_counter += 1

    finally:
        pool.close()
        pool.join()


def _init_worker(ogs: list[Orthogroup], genome_map: dict[str, Genome]):
    """
    Initializes each worker by mapping the shared memory segment.
    This runs once per CPU core.
    """
    global _WORKER_OGS, _WORKER_GENOME_MAP
    _WORKER_OGS = ogs
    _WORKER_GENOME_MAP = genome_map


def _process_og_task(og_idx: int, args) -> tuple[str, list]:
    """
    The actual computation worker.
    'og' is passed from the main process, but it references
    the static genomes now available in _SHARED_GENOME_DATA.
    """
    og = _WORKER_OGS[og_idx]

    results = og.get_refined_sogs(
        window_size=args.window, ratio_threshold=args.ratio_threshold, skip_single_ortholog=args.skip_single_orthologs
    )
    return og.og_id, results


def _calculate_optimal_chunksize(iterable_size, pool_size):
    """
    Standard heuristic used by many libraries:
    Divide the work into 4 chunks per worker to balance
    overhead vs. load balancing.
    """
    if iterable_size == 0:
        return 1

    chunksize, extra = divmod(iterable_size, pool_size * 10)
    if extra:
        chunksize += 1
    return max(1, min(chunksize, 500))


if __name__ == "__main__":
    main()
