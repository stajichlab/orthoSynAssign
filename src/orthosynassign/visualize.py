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
from pathlib import Path
from typing import TypeVar, cast

from pygenomeviz import GenomeViz
from pygenomeviz.utils import ColorCycler

from . import AUTHOR, VERSION
from ._utils import CustomHelpFormatter, VisualizeArgs, setup_logging, validate_annotations, validate_orthogroup
from .lib import Gene, get_visualize_engine, read_og_table

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
        sogs = read_og_table(sog_file, genomes)
        # Overwrite the Gene.og attribute with original og ids
        old_ogs = read_og_table(og_file, genomes)
        sogs_mapper = {og.id: idx for idx, og in enumerate(sogs)}

        # Create output directory
        if args.output is None:
            output_dir = Path("visualize_" + sog_file.stem)
        else:
            output_dir = Path(args.output)
        output_dir.mkdir(parents=True, exist_ok=True)
        logger.debug("Output directory: %s", output_dir)

        ColorCycler.set_cmap("tab20b")
        # Generate figures for each target_sog
        aligned_data_list: list[tuple[str, list[tuple[Gene, list[Gene]]]]] = []
        for target_sog in args.sog:
            if target_sog not in sogs_mapper:
                logger.warning("SOG %s not found in orthogroups.", target_sog)
                continue

            logger.debug("Collect info for %s...", target_sog)

            engine = get_visualize_engine(genomes, old_ogs, sogs)

            aligned_data_idx = engine.get_aligned_og(sogs_mapper[target_sog], args.window, args.keep_all_genes)

            aligned_data: list[tuple[Gene, list[Gene]]] = [
                (
                    genomes[genome_idx][focal_gene_idx],
                    [genomes[genome_idx][gene_idx] if gene_idx else None for gene_idx in genes],
                )
                for (genome_idx, focal_gene_idx), genes in aligned_data_idx
            ]
            aligned_data_list.append((target_sog, aligned_data))

        if args.combine and len(aligned_data_list) > 1:
            concat_aligned_data = []
            for idx, (target_sog, aligned_data) in enumerate(aligned_data_list):
                concat_aligned_data.append((target_sog, ""))
                for focal_gene, genes in aligned_data:
                    concat_aligned_data.append((focal_gene, genes))

            aligned_data_list = [("combined", concat_aligned_data)]  # Use "combined" to spy the target_sog

        for target_sog, aligned_data in aligned_data_list:
            palette = _get_palette(aligned_data)

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
        "--sog", type=str, required=True, nargs="+", help="Plot the SOG(s) refined from the previous orthosynassign analysis"
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
    opt_args.add_argument("-f", "--fmt", choices=["png", "jpg", "svg", "pdf"], default="png", help="Output image format")
    opt_args.add_argument(
        "-k", "--keep_all_genes", action="store_true", help="Keep genes that are not assigned to any orthogroup"
    )
    opt_args.add_argument(
        "-c", "--combine", action="store_true", help="Combine all SOGs in a single figure file 'combined.[fmt]'"
    )
    opt_args.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    opt_args.add_argument("-V", "--version", action="version", version=VERSION)
    opt_args.add_argument("-h", "--help", action="help", help="show this help message and exit")

    return cast(VisualizeArgs, parser.parse_args(argv))


def _get_palette(aligned_data: list[tuple[Gene, list[Gene]]]) -> dict[str, str]:
    """Prepares color palette for a SOG.

    Args:
        aligned_data (dict[Gene, list[Gene]]): A tuple containing a dictionary of aligned data of a SOG.

    Returns:
        dict[str, str]: A color palette for this SOG.
    """
    # Count OG occurrences to determine coloring
    og_counter = Counter()
    for _, win in aligned_data:
        for gene in win:
            if gene and gene.og:
                og_counter[gene.og.id] += 1

    # Build palette
    ColorCycler.reset_cycle()
    palette = {og_id: ColorCycler() for og_id, count in og_counter.items() if count > 1}

    return palette


def _render_sog_figure(aligned_windows: list[tuple[Gene, list[Gene]]], palette: dict[str, str], output_path: Path) -> None:
    """Renders the GenomeViz object and saves the file.

    Args:
        aligned_windows (dict[Gene, list[Gene]]): A dictionary of genes and their corresponding aligned windows.
        palette (dict[str, str]): A color palette for the SOG.
        output_path (Path): The path where the rendered figure should be saved.

    Returns:
        None
    """
    # Get length from the longest track
    length = max([len(gene_list) for _, gene_list in aligned_windows])

    gv = GenomeViz(fig_width=length * 1.8, track_align_type="left")

    for focal_gene, gene_list in aligned_windows:
        if isinstance(focal_gene, Gene):
            track_title = f"{focal_gene.id}\n{focal_gene.genome.name}" if isinstance(focal_gene, Gene) else focal_gene
            track = gv.add_feature_track(track_title, length * 5 - 2)
        else:
            track_title = focal_gene
            track = gv.add_feature_track(track_title, 1)

        for idx, gene in enumerate(gene_list):
            # Default values
            og_label = "None"
            text_weight = "normal"

            if gene:
                og_label = getattr(getattr(gene, "og", None), "id", "None")
                if og_label == "None":
                    fc = "#dddddd"  # Default for genes without OGs
                else:
                    fc = palette.get(og_label, "#888888")  # OGs that are not contribute to synteny

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
