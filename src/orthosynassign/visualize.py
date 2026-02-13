#!/usr/bin/env python3
"""
A companion script of orthosynassign to help visualize the refined SOGs
"""

from __future__ import annotations

import argparse
import logging
import sys
import textwrap
import traceback
from itertools import combinations
from pathlib import Path

import numpy as np
from pygenomeviz import GenomeViz
from pygenomeviz.utils import ColorCycler

from . import AUTHOR, VERSION
from ._utils import CustomHelpFormatter, setup_logging, validate_annotations, validate_orthogroup
from .lib._orthogroup import align_sog_dict, get_shared_ogs
from .lib.parsers import read_orthofinder_table

_EPILOG = textwrap.dedent(f"""\
Examples:

# Specify all bed files in a directory and plot a single SOG:
orthosynassign-vis -i orthogroup.tsv --bed *.bed --sog SOG000001.OG0000001 SOG000002.OG0000007

# Plot multiple SOGs in a single call:
orthosynassign-vis -i orthogroup.tsv --bed *.bed --sog SOG000001.OG0000001 SOG000002.OG0000007

# Specify output file name for figures:
orthosynassign -i orthogroup.tsv --bed *.bed --sog SOG000001.OG0000001 -o fig

# Specify window size applied in the previous orthosynassign analysis:
orthosynassign -i orthogroup.tsv --bed *.bed --sog SOG000001.OG0000001 -o fig -w 10

# With verbose output:
orthosynassign-vis -i orthogroup.tsv --bed *.bed --sog SOG000001.OG0000001 -v

Written by {AUTHOR}
""")


def main(args: list[str] | None = None) -> int:
    """OrthoSynAssign visualization script."""
    # Parse arguments
    args = _parse_arguments() or sys.argv[1:]

    # Setup logging
    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)

    logger.info("Starting Visualize")
    logger.debug(f"Command: {' '.join(sys.argv)}")

    try:
        # Validate inputs
        annotations = validate_annotations(args, logger)
        og_file = validate_orthogroup(args)

        # Read gff
        genomes = {}
        for annotation in annotations:
            genome = annotation.parse()
            genomes[genome.sample_name] = genome

        # Read orthogroups
        logger.info(f"Reading orthogroup data from: {og_file}")
        orthogroups = read_orthofinder_table(og_file, genomes)
        logger.debug(f"Loaded {len(orthogroups)} orthogroups")

        # Create output directory
        if args.output is None:
            parent_dir = og_file.parent
            output_dir = parent_dir / ("visualize_" + og_file.stem)
        else:
            output_dir = Path(args.output)
        output_dir.mkdir(parents=True, exist_ok=True)
        logger.debug(f"Output directory: {output_dir}")

        # Perform synteny analysis
        og_id_list = [og.og_id for og in orthogroups]
        ColorCycler.set_cmap("tab10")
        for target_sog in args.sog:
            try:
                idx = og_id_list.index(target_sog)
            except ValueError:
                logger.warning("SOG %s not found in orthogroups.", target_sog)
                continue

            logger.debug("Generate figure for %s...", target_sog)
            og_windows = {gene: [] for gene in orthogroups[idx]}
            for gene_A, gene_B in combinations(orthogroups[idx], 2):
                og_windows[gene_A] = _merge_by_overlap(
                    og_windows[gene_A],
                    gene_A.genome.get_window(
                        gene_A, get_shared_ogs(gene_A.genome, gene_B.genome), window_size=args.window, mode="validate"
                    ),
                )
                og_windows[gene_B] = _merge_by_overlap(
                    og_windows[gene_B],
                    gene_B.genome.get_window(
                        gene_B, get_shared_ogs(gene_A.genome, gene_B.genome), window_size=args.window, mode="validate"
                    ),
                )

            aligned_og_windows = align_sog_dict(og_windows)
            uniq_ogs = {}
            for win in list(aligned_og_windows.values()):
                for gene in win:
                    if gene is not None:
                        if gene.orthogroup:
                            uniq_ogs[(gene.orthogroup.og_id.split(".")[-1])] = ""

            ColorCycler.reset_cycle()

            for idx, og in enumerate(uniq_ogs):
                uniq_ogs[og] = ColorCycler()

            length = len(list(aligned_og_windows.values())[0])
            gv = GenomeViz(fig_width=length * 1.8, track_align_type="center")
            for focal_gene, gene_list in aligned_og_windows.items():
                track_title = f"{focal_gene.id}\n{focal_gene.genome.sample_name}"
                track = gv.add_feature_track(track_title, length * 5 - 2)

                for idx, gene in enumerate(gene_list):
                    text_weight = "normal"
                    if gene is not None:
                        if gene.orthogroup:
                            if gene == focal_gene:
                                text_weight = "bold"
                            og = gene.orthogroup.og_id.split(".")[-1]
                            fc = uniq_ogs[og]
                        else:
                            og = "None"
                            fc = "#777777"
                        start = idx * 5
                        track.add_feature(
                            start,
                            start + 3,
                            plotstyle="bigrbox",
                            label="\n".join(textwrap.wrap(gene.id, width=round(length * 0.8))),
                            fc=fc,
                            text_kws={"weight": text_weight, "rotation": 0, "hpos": "center"},
                        )
                        track.add_text(start + 1.5, og, vpos="bottom", hpos="center", rotation=0, weight=text_weight)

            output_file_path = output_dir / f"{target_sog}.png"
            gv.savefig(output_file_path)
            logger.info(f"Figure saved to {output_file_path}")

        logger.info("Visualize completed successfully")

    except KeyboardInterrupt:
        logger.warning("Terminated by user.")
        return 1

    except Exception as e:
        logger.error(f"An error occurred: {e}")
        logger.debug("%s", traceback.format_exc())
        return 1

    return 0


def _parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=CustomHelpFormatter,
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

    req_args.add_argument(
        "--bed",
        type=Path,
        required=True,
        metavar=("file", "files"),
        nargs="+",
        help="Path of BED formatted genome annotation files",
    )

    req_args.add_argument(
        "--sog", type=str, required=True, nargs="+", help="Plot the SOG of the previous orthosynassign analysis"
    )

    opt_args = parser.add_argument_group("Options")
    opt_args.add_argument(
        "-o",
        "--output",
        dest="output",
        type=Path,
        help="Output directory (default: visualize_[og_input])",
    )

    opt_args.add_argument(
        "-w",
        "--window",
        dest="window",
        type=int,
        default=8,
        help="The window size applied to the previous orthosynassign analysis",
    )
    opt_args.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    opt_args.add_argument("-V", "--version", action="version", version=VERSION)
    opt_args.add_argument("-h", "--help", action="help", help="show this help message and exit")

    return parser.parse_args()


def _merge_by_overlap(list_a, list_b):
    """
    Merges list_b into list_a by finding the first overlap index.
    """
    arr_a = np.array(list_a)
    arr_b = np.array(list_b)

    if arr_a.size == 0:
        return list_b

    # Find indices in A where the first element of B appears
    matches = np.where(arr_a == arr_b[0])[0]

    if matches.size > 0:
        # We take the first match found
        overlap_start_idx = matches[0]
        # Concatenate: everything from A up to the overlap, plus all of B
        return list(np.concatenate([arr_a[:overlap_start_idx], arr_b]))

    # No overlap found, just join them (or handle as a gap)
    return list(np.concatenate([arr_a, arr_b]))
