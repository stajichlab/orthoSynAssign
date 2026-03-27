#!/usr/bin/env python3
"""
A companion script of orthosynassign to help visualize the refined SOGs
"""

from __future__ import annotations

import argparse
import logging
import sys
import textwrap
from collections import Counter
from itertools import combinations
from pathlib import Path
from typing import TYPE_CHECKING, TypeVar, cast

from pygenomeviz import GenomeViz
from pygenomeviz.utils import ColorCycler

from . import AUTHOR, VERSION
from ._utils import CustomHelpFormatter, VisualizeArgs, setup_logging, validate_annotations, validate_orthogroup
from .lib import align_sog_dict, get_shared_ogs, read_og_table

if TYPE_CHECKING:
    from .lib import Gene, Orthogroup

_EPILOG = textwrap.dedent(f"""\
Examples:

# Plot multiple SOGs in a single call:
orthosynassign-vis --og_file orthogroup.tsv --sog_file Refined_SOGs.tsv --bed *.bed --sog SOG000001.OG0000001 SOG000002.OG0000007

# Specify output directory for figures:
orthosynassign --og_file orthogroup.tsv --sog_file Refined_SOGs.tsv --bed *.bed --sog SOG000001.OG0000001 -o fig

# Specify window size applied in the previous orthosynassign analysis:
orthosynassign --og_file orthogroup.tsv --sog_file Refined_SOGs.tsv --bed *.bed --sog SOG000001.OG0000001 -w 10

# Save in svg format:
orthosynassign-vis --og_file orthogroup.tsv --sog_file Refined_SOGs.tsv --bed *.bed --sog SOG000001.OG0000001 -f svg

# With verbose output:
orthosynassign-vis --og_file orthogroup.tsv --sog_file Refined_SOGs.tsv --bed *.bed --sog SOG000001.OG0000001 -v

Written by {AUTHOR}
""")

_T = TypeVar("Type")


def run_cli() -> None:
    """Runs the orthosynassign-vis CLI entry point."""
    parsed: VisualizeArgs = _parse_arguments(sys.argv[1:])
    sys.exit(main(parsed))


def main(args: VisualizeArgs) -> int:
    """Main entry point for orthosynassign-vis.

    Args:
        args (VisualizeArgs): Parsed command line arguments.

    Returns:
        int: Exit code.
    """
    # Setup logging
    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)

    logger.info("Starting Visualize")
    logger.debug("Command: %s ".join(sys.argv))

    try:
        # Validate inputs
        annotations = validate_annotations(args)
        og_file = validate_orthogroup(args.og_file)
        sog_file = validate_orthogroup(args.sog_file)

        # Read gff
        genomes = []
        for annotation in annotations:
            genome = annotation.parse()
            genomes.append(genome)

        # Read orthogroups
        logger.info("Reading orthogroup data from: %s", og_file)
        refined_orthogroups = {og.id: og for og in read_og_table(sog_file, genomes)}
        # Overwrite the Gene.og attribute with original og ids
        _ = read_og_table(og_file, genomes)

        # Create output directory
        if args.output is None:
            output_dir = Path("visualize_" + sog_file.stem)
        else:
            output_dir = Path(args.output)
        output_dir.mkdir(parents=True, exist_ok=True)
        logger.debug("Output directory: %s", output_dir)

        ColorCycler.set_cmap("tab20b")
        # Generate figures for each target_sog
        for target_sog in args.sog:
            if target_sog not in refined_orthogroups:
                logger.warning("SOG %s not found in orthogroups.", target_sog)
                continue

            logger.debug("Generate figure for %s...", target_sog)

            aligned_data, palette = _prepare_sog_visualization_data(
                refined_orthogroups[target_sog], window_size=args.window, keep_all_genes=args.keep_all_genes
            )

            output_file_path = output_dir / f"{target_sog}.{args.fmt}"
            _render_sog_figure(aligned_data, palette, output_file_path)

            logger.info("Figure saved to %s", output_file_path)

        logger.info("Visualize completed successfully")

    except KeyboardInterrupt:
        logger.warning("Terminated by user.")
        return 130

    except FileNotFoundError as e:
        logger.error("An error occurred: %s", e)
        logger.debug("Traceback details:", exc_info=True)
        return 2

    except Exception as e:
        logger.error("An error occurred: %s", e)
        logger.debug("Traceback details:", exc_info=True)
        return 1

    return 0


def _parse_arguments(argv=None) -> VisualizeArgs:
    """Parse command line arguments.

    Args:
        argv (list of str, optional): The list of arguments to parse. Defaults to sys.argv[1:].

    Returns:
        VisualizeArgs: Parsed command line arguments.
    """
    parser = argparse.ArgumentParser(
        description=__doc__,
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
        help="Path to the original orthogroups.tsv file",
    )
    req_args.add_argument(
        "--sog_file",
        type=Path,
        required=True,
        help="Path to the refined orthogroups.tsv file",
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
        "-w",
        "--window",
        type=int,
        default=8,
        help="The window size applied to the previous orthosynassign analysis",
    )
    opt_args.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Output directory (default: visualize_[sog_file])",
    )
    opt_args.add_argument("-f", "--fmt", choices=["png", "jpg", "svg", "pdf"], default="png", help="Output image format.")
    opt_args.add_argument(
        "-k", "--keep_all_genes", action="store_true", help="Keep genes that are not assigned to any orthogroup"
    )
    opt_args.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    opt_args.add_argument("-V", "--version", action="version", version=VERSION)
    opt_args.add_argument("-h", "--help", action="help", help="show this help message and exit")

    return cast(VisualizeArgs, parser.parse_args(argv))


def _prepare_sog_visualization_data(
    og: Orthogroup, *, window_size: int, keep_all_genes=False
) -> tuple[dict[Gene, list[Gene]], dict[str, str]]:
    """Prepares aligned windows and a color palette for a SOG.

    Args:
        og (Orthogroup): The orthogroup to prepare visualization data for.
        window_size (int): The size of the window around each gene in the SOG.
        keep_all_genes (bool, optional): Whether to include all genes or only those assigned to an orthogroup. Defaults to False.

    Returns:
        tuple[dict[Gene, list[Gene]], dict[str, str]]: A tuple containing a dictionary of aligned windows and a color palette for
        the SOG.
    """
    og_windows: dict[Gene, list[Gene]] = _get_genes_from_window(og, window=window_size, keep_all_genes=keep_all_genes)
    aligned_windows: dict[Gene, list[Gene]] = align_sog_dict(og_windows)

    # Count OG occurrences to determine coloring
    og_counter = Counter()
    for win in aligned_windows.values():
        for gene in win:
            if gene and gene.og:
                og_counter[gene.og.id] += 1

    # Build palette
    ColorCycler.reset_cycle()
    palette = {og_id: ColorCycler() for og_id, count in og_counter.items() if count > 1}

    return aligned_windows, palette


def _get_genes_from_window(og: Orthogroup, *, window: int, keep_all_genes=False) -> dict[Gene, list[Gene]]:
    """Retrieves a dictionary of genes and their corresponding windows from the given orthogroup.

    Args:
        og (Orthogroup): The orthogroup to retrieve genes from.
        window (int): The size of the window around each gene in the SOG.
        keep_all_genes (bool, optional): Whether to include all genes or only those assigned to an orthogroup. Defaults to False.

    Returns:
        dict[Gene, list[Gene]]: A dictionary where the keys are genes and the values are lists of genes within their respective
        windows.
    """
    og_windows_idx: dict[Gene, list[int]] = {gene: [gene.index, gene.index] for gene in og if gene.index}
    for gene_A, gene_B in combinations(og, 2):
        if not gene_A.genome or not gene_B.genome:
            continue
        win_A = gene_A.genome.get_window(gene_A, get_shared_ogs(gene_A.genome, gene_B.genome), window_size=window)
        win_B = gene_B.genome.get_window(gene_B, get_shared_ogs(gene_A.genome, gene_B.genome), window_size=window)
        for gene, win in zip((gene_A, gene_B), (win_A, win_B)):
            start, end = og_windows_idx[gene]
            new_start, new_end = win[0].index, win[-1].index
            start = new_start if new_start < start else start
            end = new_end if new_end > end else end
            og_windows_idx[gene] = [start, end]

    og_windows: dict[Gene, list[Gene]] = {gene: [] for gene in og}
    for gene, (start, end) in og_windows_idx.items():
        continuous_slice: list[Gene] = gene.genome[start : end + 1]
        if keep_all_genes:
            og_windows[gene] = list(continuous_slice)
        else:
            # Filter: Keep if gene has an OG or is the focal gene itself
            og_windows[gene] = [g for g in continuous_slice if g.og is not None or g == gene.representative]

    return og_windows


def _render_sog_figure(aligned_windows: dict[Gene, list[Gene]], palette: dict[str, str], output_path: Path) -> None:
    """Renders the GenomeViz object and saves the file.

    Args:
        aligned_windows (dict[Gene, list[Gene]]): A dictionary of genes and their corresponding aligned windows.
        palette (dict[str, str]): A color palette for the SOG.
        output_path (Path): The path where the rendered figure should be saved.

    Returns:
        None
    """
    # Get length from the first available track
    first_track = next(iter(aligned_windows.values()))
    length = len(first_track)

    gv = GenomeViz(fig_width=length * 1.8, track_align_type="center")

    for focal_gene, gene_list in aligned_windows.items():
        track_title = f"{focal_gene.id}\n{focal_gene.genome.name}"
        track = gv.add_feature_track(track_title, length * 5 - 2)

        for idx, gene in enumerate(gene_list):
            # Default values
            og_label = "None"
            fc = "#cccccc"  # Default for genes without OGs
            text_weight = "normal"

            if gene:
                og_label = getattr(getattr(gene, "og", None), "id", "None")
                fc = palette.get(og_label, "#777777")

                if gene == focal_gene.representative:
                    text_weight = "bold"
                    fc = "#fcfc42"  # Highlight focal gene

                start = idx * 5
                track.add_feature(
                    start,
                    start + 3,
                    plotstyle="bigrbox",
                    label="\n".join(textwrap.wrap(gene.id, width=12)),
                    fc=fc,
                    text_kws={"weight": text_weight, "rotation": 0, "hpos": "center"},
                )
                track.add_text(start + 1.5, og_label, vpos="bottom", hpos="center", rotation=0, weight=text_weight)

    gv.savefig(output_path)


if __name__ == "__main__":
    run_cli()
