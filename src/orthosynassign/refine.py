#!/usr/bin/env python3
"""
The orthosynassign CLI entry point.
"""

from __future__ import annotations

import argparse
import logging
import multiprocessing
import os
import sys
import textwrap
import time
from functools import partial
from pathlib import Path
from typing import TYPE_CHECKING, Iterator, cast

from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm

from . import AUTHOR, VERSION
from . import __doc__ as _module_doc
from ._utils import CustomHelpFormatter, RefineArgs, setup_logging, validate_annotations, validate_orthogroup
from .lib.parsers import read_orthogroup_table, save_results_tsv

if TYPE_CHECKING:
    from .lib._gene import Gene, Genome
    from .lib._orthogroup import Orthogroup

_EPILOG = textwrap.dedent(f"""\
Examples:

# Specify bed files separately:
orthosynassign --og_file orthogroup.tsv --bed file1.bed file2.bed file3.bed

# Specify all bed files in a directory and processed in parallel with 6 CPUs:
orthosynassign --og_file orthogroup.tsv --bed *.bed -t 6

# Specify output file name for results:
orthosynassign --og_file orthogroup.tsv --bed *.bed -o Refined_SOGs.tsv

# Specify window size and ratio threshold:
orthosynassign --og_file orthogroup.tsv --bed *.bed -w 10 -r 0.8

# With verbose output:
orthosynassign --og_file orthogroup.tsv --bed *.bed -v

Written by {AUTHOR}
""")

_WORKER_OGS: list[Orthogroup] = []
_WORKER_GENOME_MAP: dict[str, Genome] = {}
_AVAIL_CPUS = int(os.environ.get("SLURM_CPUS_ON_NODE", os.cpu_count()))


def run_cli() -> None:
    """Runs the orthoSynAssign CLI entry point."""
    parsed: RefineArgs = _parse_arguments(sys.argv[1:])
    sys.exit(main(parsed))


def main(args: RefineArgs) -> int:
    """Main entry point for orthoSynAssign.

    Args:
        args (RefineArgs): Parsed command line arguments.

    Returns:
        int: Exit code.
    """
    # Setup logging
    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)

    logger.info("Starting orthoSynAssign")
    logger.debug("Command: %s", " ".join(sys.argv))

    output_path = Path(args.output)
    tmp_output = output_path.with_suffix(output_path.suffix + ".tmp")

    try:
        # Validate inputs
        annotations = validate_annotations(args)
        og_file = validate_orthogroup(args.og_file)

        # Create output directory
        output_dir = output_path.parent
        output_dir.mkdir(parents=True, exist_ok=True)
        logger.info("Creating output directory: %s", output_dir)

        # Read gff
        genomes = {}
        for annotation in annotations:
            genome = annotation.parse()
            genomes[genome.name] = genome

        # Read orthogroup
        logger.info("Reading orthogroup data from: %s", og_file)
        orthogroups = read_orthogroup_table(og_file, genomes)

        # Perform synteny analysis
        logger.info("Refining orthogroups by pairwise synteny analysis.")
        results_stream = _generate_sog_results(orthogroups, genomes, args, cpus=args.threads)

        save_results_tsv(results_stream, list(genomes.keys()), tmp_output)
        tmp_output.replace(output_path)
        logger.info("Refinement complete. Results saved to %s", args.output)

        logger.info("orthoSynAssign completed successfully")

    except KeyboardInterrupt:
        logger.warning("Terminated by user.")
        return 1

    except Exception as e:
        logger.error("An error occurred: %s", e)
        logger.debug("Traceback details:", exc_info=True)
        return 1

    finally:
        if tmp_output.exists():
            tmp_output.unlink()

    return 0


def _parse_arguments(argv=None) -> RefineArgs:
    """Parse command line arguments.

    Args:
        argv (list of str, optional): The list of arguments to parse. Defaults to sys.argv[1:].

    Returns:
        RefineArgs: Parsed command line arguments.
    """
    parser = argparse.ArgumentParser(
        description=_module_doc,
        formatter_class=CustomHelpFormatter,
        epilog=_EPILOG,
        add_help=False,
    )
    req_args = parser.add_argument_group("Required arguments")
    # OrthoFinder input
    req_args.add_argument(
        "--og_file",
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
        type=int,
        default=8,
        help="Controls how many total genes are considered when determining synteny for a single gene",
    )

    opt_args.add_argument(
        "-r",
        "--ratio_threshold",
        dest="threshold",
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
        type=Path,
        default="Refined_SOGs-%s.tsv" % time.strftime("%Y%m%d-%H%M%S", time.gmtime()),
        help="Output of results (default: Refined_SOGs-[YYYYMMDD-HHMMSS].tsv (UTC timestamp))",
    )
    opt_args.add_argument("-t", "--threads", type=int, default=min(_AVAIL_CPUS, 4), help="Number of cpus to use")
    opt_args.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    opt_args.add_argument("-V", "--version", action="version", version=VERSION)
    opt_args.add_argument("-h", "--help", action="help", help="show this help message and exit")

    return cast(RefineArgs, parser.parse_args(argv))


def _generate_sog_results(
    orthogroups: list[Orthogroup], genome_data: dict[str, Genome], args: RefineArgs, *, cpus: int = 1
) -> Iterator[tuple[str, dict[Genome, list[Gene]]]]:
    """
    Processes orthogroups and yields results one by one.

    Args:
        orthogroups (list[Orthogroup]): The list of orthogroups to process.
        genome_data (dict[str, Genome]): A dictionary mapping genome names to their respective Genome objects.
        args (RefineArgs): Command-line arguments.
        cpus (int, optional): Number of CPUs to use for parallel processing. Defaults to 1.

    Yields:
        Iterator[tuple[str, dict[Genome, list[Gene]]]]: An iterator yielding tuples of SOG ID and a dictionary mapping genomes to
        lists of genes.
    """
    global_sog_counter = 1
    total_ogs = len(orthogroups)
    indices = list(range(total_ogs))
    desc = f"Refining with {cpus} cpu{'s' if cpus > 1 else ''}"

    tqdm_kwargs = {"desc": desc, "unit": "og", "total": total_ogs, "mininterval": 1, "ascii": " #"}

    with logging_redirect_tqdm():
        if cpus == 1:
            _init_worker(orthogroups, genome_data)

            for index in tqdm(indices, **tqdm_kwargs):
                _, refined_sogs = _process_og_task(index, args)

                for sog_dict in refined_sogs:
                    sog_id = f"SOG{global_sog_counter:06d}.{orthogroups[index].id}"
                    yield sog_id, sog_dict
                    global_sog_counter += 1
            return

        # Parallel logic
        ctx = multiprocessing.get_context("spawn")
        opt_chunksize = _calculate_optimal_chunksize(total_ogs, cpus)
        pool = ctx.Pool(processes=cpus, initializer=_init_worker, initargs=(orthogroups, genome_data))

        try:
            worker_func = partial(_process_og_task, args=args)

            # imap returns an iterator, tqdm wraps it
            pbar = tqdm(pool.imap(worker_func, indices, chunksize=opt_chunksize), **tqdm_kwargs)

            for og_id, refined_sogs in pbar:
                for sog_dict in refined_sogs:
                    sog_id_str = f"SOG{global_sog_counter:06d}.{og_id}"
                    yield sog_id_str, sog_dict
                    global_sog_counter += 1

        finally:
            pool.close()
            pool.join()


def _init_worker(ogs: list[Orthogroup], genome_map: dict[str, Genome]) -> None:
    """Initializes each worker by setting up the global variables that will be used for processing orthogroups. This runs once per
    CPU core."""
    global _WORKER_OGS, _WORKER_GENOME_MAP
    _WORKER_OGS = ogs
    _WORKER_GENOME_MAP = genome_map


def _process_og_task(og_idx: int, args: RefineArgs) -> tuple[str, list[dict[Genome, list[Gene]]]]:
    """The actual computation worker that processes a single orthogroup.

    'og' is passed from the main process, but it references the static genomes now available in _WORKER_GENOME_MAP.

    Args:
        og_idx (int): Index of the orthogroup to process.
        args (RefineArgs): Command-line arguments for processing.

    Returns:
        tuple[str, list[dict[Genome, list[Gene]]]]: A tuple containing
            the orthogroup ID and a list of dictionaries mapping genomes
            to lists of refined genes.
    """
    og = _WORKER_OGS[og_idx]

    results = og.get_refined_sogs(window_size=args.window, ratio_threshold=args.threshold)
    return og.id, results


def _calculate_optimal_chunksize(iterable_size: int, pool_size: int) -> int:
    """
    Calculate the optimal chunk size for dividing work among workers.

    Args:
        iterable_size (int): The total number of items to be processed.
        pool_size (int): The number of worker processes available.

    Returns:
        int: The calculated optimal chunk size. This is intended to balance the overhead of task distribution with the load
            balancing across workers. A standard heuristic used by many libraries suggests dividing the work into 4 chunks per
            worker.
    """
    if iterable_size == 0:
        return 1

    chunksize, extra = divmod(iterable_size, pool_size * 10)
    if extra:
        chunksize += 1
    return max(1, min(chunksize, 1000))


if __name__ == "__main__":
    run_cli()
