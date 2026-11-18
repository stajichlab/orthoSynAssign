#!/usr/bin/env python3
"""
A companion script of orthosynassign to help visualize the refined SOGs
"""

from __future__ import annotations

import argparse
import logging
import sys
import textwrap
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from itertools import chain, combinations
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

# Collapse tracks whose windows overlap instead of only tandem copies:
orthosynassign-vis --og_file orthogroup.tsv --sog_file Refined_SOGs.tsv --bed *.bed --sog SOG000001.OG0000001 --collapse window

# With verbose output:
orthosynassign-vis --og_file orthogroup.tsv --sog_file Refined_SOGs.tsv --bed *.bed --sog SOG000001.OG0000001 -v

Written by {AUTHOR}
""")

_T = TypeVar("Type")

#: Color used to highlight the focal gene(s) of a track.
_FOCAL_COLOR = "#fcfc42"

#: Maximum number of characters per line of a gene label.
_LABEL_WIDTH = 12

#: Shortest last line a wrapped gene label may end with.
_LABEL_MIN_TAIL = 2

#: Color of the bar that marks the end of a contig.
_CONTIG_END_COLOR = "#444444"

#: Available modes for --collapse.
COLLAPSE_MODES = ("tandem", "window", "none")


@dataclass(frozen=True)
class Track:
    """One row of a figure, fully shaped and ready to draw.

    Every decision taken from the command line options (`-k/--keep_all_genes`, `--collapse`)
    is resolved while the tracks are built, so rendering only has to walk the genes.

    Attributes:
        title (str): The text drawn to the left of the track.
        genes (list[Gene | None]): The genes of the aligned window, `None` where the track has
            no gene in that column.
        focal_ids (frozenset[str]): The IDs of the focal genes the track stands for, including
            the genes they represent. These are highlighted with the focal gene color.
        at_contig_start (bool): The first gene of the track is the first gene of its contig, so
            the window is cut short by the contig and not by the window size.
        at_contig_end (bool): The last gene of the track is the last gene of its contig.
    """

    title: str
    genes: list[Gene | None] = field(default_factory=list)
    focal_ids: frozenset[str] = frozenset()
    at_contig_start: bool = False
    at_contig_end: bool = False

    @property
    def is_separator(self) -> bool:
        """Reports whether the track is an empty spacer between two SOGs.

        Returns:
            bool: True if the track holds no gene.
        """
        return not self.genes


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
        tracks_per_sog: list[tuple[str, list[Track]]] = []
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
            tracks = _build_tracks(aligned_data, args.collapse)
            if len(tracks) < len(aligned_data):
                logger.info(
                    "%s: collapsed %d focal genes into %d tracks (--collapse %s)",
                    target_sog,
                    len(aligned_data),
                    len(tracks),
                    args.collapse,
                )
            tracks_per_sog.append((target_sog, tracks))

        if args.combine and len(tracks_per_sog) > 1:
            concat_tracks: list[Track] = []
            for target_sog, tracks in tracks_per_sog:
                # A track without genes acts as the separator between two SOGs
                concat_tracks.append(Track(title=target_sog))
                concat_tracks.extend(tracks)

            tracks_per_sog = [("combined", concat_tracks)]  # Use "combined" to spy the target_sog

        for target_sog, tracks in tracks_per_sog:
            palette = _get_palette(tracks)

            output_file_path = output_dir / f"{target_sog}.{args.fmt}"
            _render_sog_figure(tracks, palette, output_file_path)

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
        "--collapse",
        choices=list(COLLAPSE_MODES),
        default="tandem",
        help=(
            "Collapse focal genes of the same genome into a single track: "
            "'tandem' collapses focal genes that sit next to each other, "
            "'window' collapses focal genes whose windows overlap, "
            "'none' disables collapsing"
        ),
    )
    opt_args.add_argument(
        "-c", "--combine", action="store_true", help="Combine all SOGs in a single figure file 'combined.[fmt]'"
    )
    opt_args.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    opt_args.add_argument("-V", "--version", action="version", version=VERSION)
    opt_args.add_argument("-h", "--help", action="help", help="show this help message and exit")

    return cast(VisualizeArgs, parser.parse_args(argv))


def _build_tracks(aligned_data: list[tuple[Gene, list[Gene]]], mode: str = "tandem") -> list[Track]:
    """Shapes the aligned windows of a SOG into the tracks to draw.

    Args:
        aligned_data (list[tuple[Gene, list[Gene]]]): Focal genes and their aligned windows.
        mode (str): The `--collapse` mode.

    Returns:
        list[Track]: One track per figure row.
    """
    rows = _collapse_tracks(aligned_data, mode)
    edges = [_contig_edges(window) for _, window, _ in rows]

    # The marks are drawn beside the first and last gene, which needs a column to live in
    left_margin = [None] if any(at_start for at_start, _ in edges) else []
    right_margin = [None] if any(at_end for _, at_end in edges) else []

    tracks = []
    for (focal_gene, window, members), (at_start, at_end) in zip(rows, edges):
        # A member may stand for tandem in-paralogs of its own, which are focal genes as well
        focal_ids = {gene.id for gene in members} | {gene.representative.id for gene in members}
        tracks.append(
            Track(
                title=f"{_track_label(members)}\n{focal_gene.genome.name}",
                genes=left_margin + window + right_margin,
                focal_ids=frozenset(focal_ids),
                at_contig_start=at_start,
                at_contig_end=at_end,
            )
        )

    return tracks


def _contig_edges(window: list[Gene]) -> tuple[bool, bool]:
    """Reports whether a window reaches the start or the end of its contig.

    A window stops either because the window size is used up or because the contig ends, and
    the figure would look the same in both cases. The engine does not report which of the two
    happened, but the genes of a contig are consecutive in the genome, so the neighbours of
    the first and the last gene of the window answer it. A circular genome has no end.

    Args:
        window (list[Gene]): The genes of the window, padding included.

    Returns:
        tuple[bool, bool]: Whether the window starts at the first gene of its contig and
        whether it stops at the last one.
    """
    genes = [gene for gene in window if gene]
    if not genes or genes[0].genome is None or genes[0].genome.is_circular:
        return (False, False)

    return (_at_contig_edge(genes[0], -1), _at_contig_edge(genes[-1], 1))


def _at_contig_edge(gene: Gene, step: int) -> bool:
    """Reports whether a gene is the first or the last gene of its contig.

    Args:
        gene (Gene): The gene to test.
        step (int): -1 to test the start of the contig, 1 to test its end.

    Returns:
        bool: True if the genome holds no gene of the same contig on that side.
    """
    genome = gene.genome
    neighbour = gene.index + step
    if neighbour < 0 or neighbour >= len(genome):
        return True

    return genome[neighbour].seqid != gene.seqid


def _collapse_tracks(
    aligned_data: list[tuple[Gene, list[Gene]]], mode: str = "tandem"
) -> list[tuple[Gene, list[Gene], list[Gene]]]:
    """Collapses focal genes of the same genome into a single track.

    Args:
        aligned_data (list[tuple[Gene, list[Gene]]]): Focal genes and their aligned windows.
        mode (str): "tandem" collapses focal genes that are adjacent on the same contig,
            "window" collapses focal genes whose windows overlap, "none" disables collapsing.

    Returns:
        list[tuple[Gene, list[Gene], list[Gene]]]: The focal gene, the aligned window and the
        focal genes represented by each kept track.

    Raises:
        ValueError: If an unknown mode is given.
    """
    if mode not in COLLAPSE_MODES:
        raise ValueError(f"Unknown collapse mode: {mode}. Choose from {', '.join(COLLAPSE_MODES)}.")

    if mode == "none" or len(aligned_data) < 2:
        covers = [(idx, [idx]) for idx in range(len(aligned_data))]
    else:
        covers = []
        for group in _group_tracks(aligned_data, mode):
            covers.extend(_split_group(group, aligned_data))

    # Keep the tracks in the order the engine reports them, that is by genome and gene index
    collapsed = [
        (
            aligned_data[kept][0],
            _merge_windows([aligned_data[idx][1] for idx in members]),
            [aligned_data[idx][0] for idx in members],
        )
        for kept, members in sorted(covers)
    ]

    return _align_tracks(collapsed)


def _split_group(group: list[int], aligned_data: list[tuple[Gene, list[Gene]]]) -> list[tuple[int, list[int]]]:
    """Splits a group into as few tracks as possible without hiding a focal gene.

    Grouping is transitive, so a chain of overlapping windows can reach much further than any
    single window: the focal genes at one end of the chain are then missing from the window
    kept for the other end, and would not be drawn at all. The group is therefore covered
    greedily, taking the window that shows the most focal genes still to be placed, until
    every focal gene of the group is shown by the track it was merged into.

    Args:
        group (list[int]): The track indices to cover.
        aligned_data (list[tuple[Gene, list[Gene]]]): Focal genes and their aligned windows.

    Returns:
        list[tuple[int, list[int]]]: The kept track and the track indices it stands for, one
        entry per track the group is split into.
    """
    remaining = sorted(group)
    covers: list[tuple[int, list[int]]] = []

    while remaining:
        # A window always holds its own focal gene, so every round places at least one track
        kept = max(
            remaining,
            key=lambda idx: (
                len(_covered_by(idx, remaining, aligned_data)),
                sum(1 for gene in aligned_data[idx][1] if gene),
                -idx,
            ),
        )
        members = _covered_by(kept, remaining, aligned_data)
        covers.append((kept, members))
        remaining = [idx for idx in remaining if idx not in members]

    return covers


def _covered_by(kept: int, candidates: list[int], aligned_data: list[tuple[Gene, list[Gene]]]) -> list[int]:
    """Collects the tracks whose focal gene is shown by the window of another track.

    Args:
        kept (int): The index of the track that would be kept.
        candidates (list[int]): The track indices to test.
        aligned_data (list[tuple[Gene, list[Gene]]]): Focal genes and their aligned windows.

    Returns:
        list[int]: The indices of the tracks the kept window can stand for.
    """
    shown = _window_keys(aligned_data[kept][1])
    return [idx for idx in candidates if _gene_key(aligned_data[idx][0]) in shown]


def _merge_windows(windows: list[list[Gene]]) -> list[Gene]:
    """Merges the windows of the tracks that are drawn as one into a single window.

    Keeping only the window of the track that is kept would drop the genes that the other
    windows show beyond its ends. The windows of a merged track overlap and each of them is a
    run of consecutive genes, so the widest one is extended with what the others show before
    its first and after its last gene. Working from the overlap, rather than from the gene
    index, also holds for a window that wraps around the origin of a circular genome.

    Args:
        windows (list[list[Gene]]): The aligned windows to merge, padding included.

    Returns:
        list[Gene]: The genes of the merged window, in the order they are drawn.
    """
    runs = sorted(([gene for gene in window if gene] for window in windows), key=len, reverse=True)
    merged = runs[0]

    for run in runs[1:]:
        shared = [idx for idx, gene in enumerate(run) if _gene_key(gene) in {_gene_key(g) for g in merged}]
        if not shared:
            # Not reachable for a merged track, whose windows always overlap
            continue
        merged = run[: shared[0]] + merged + run[shared[-1] + 1 :]

    return merged


def _align_tracks(
    tracks: list[tuple[Gene, list[Gene], list[Gene]]],
) -> list[tuple[Gene, list[Gene], list[Gene]]]:
    """Lines the tracks up on their first focal gene by padding them with empty columns.

    The engine aligns the windows it returns, but collapsing changes both which tracks are
    drawn and how far each of them reaches, so the padding is rebuilt here. Only the columns
    needed by the tracks that are actually drawn are kept, which is what stops a removed track
    from leaving blank space behind.

    Args:
        tracks (list[tuple[Gene, list[Gene], list[Gene]]]): The tracks and their merged windows.

    Returns:
        list[tuple[Gene, list[Gene], list[Gene]]]: The tracks, their windows padded with `None`
        to a common length.
    """
    if not tracks:
        return tracks

    anchors = []
    for _, window, members in tracks:
        keys = [_gene_key(gene) for gene in window]
        # The first focal gene of a track is the one its neighbours are aligned on
        anchors.append(min((keys.index(_gene_key(gene)) for gene in members if _gene_key(gene) in keys), default=0))

    max_anchor = max(anchors)
    padded = [[None] * (max_anchor - anchor) + list(window) for (_, window, _), anchor in zip(tracks, anchors)]
    length = max(len(window) for window in padded)

    return [
        (focal_gene, window + [None] * (length - len(window)), members)
        for (focal_gene, _, members), window in zip(tracks, padded)
    ]


def _gene_key(gene: Gene) -> tuple[str, int]:
    """Identifies a gene by its position in its genome.

    Args:
        gene (Gene): The gene to identify.

    Returns:
        tuple[str, int]: The contig and index of the gene.
    """
    return (gene.seqid, gene.index)


def _group_tracks(aligned_data: list[tuple[Gene, list[Gene]]], mode: str) -> list[list[int]]:
    """Groups the indices of the tracks that should be collapsed together.

    Only focal genes of the same genome are ever grouped. Grouping is transitive, so a run of
    three tandem copies ends up in a single group.

    Args:
        aligned_data (list[tuple[Gene, list[Gene]]]): Focal genes and their aligned windows.
        mode (str): Either "tandem" or "window".

    Returns:
        list[list[int]]: Groups of track indices, each sorted, ordered by their first track.
    """
    parent = list(range(len(aligned_data)))

    def find(idx: int) -> int:
        while parent[idx] != idx:
            parent[idx] = parent[parent[idx]]
            idx = parent[idx]
        return idx

    def union(a: int, b: int) -> None:
        root_a, root_b = find(a), find(b)
        if root_a != root_b:
            parent[max(root_a, root_b)] = min(root_a, root_b)

    by_genome: dict[str, list[int]] = defaultdict(list)
    for idx, (focal_gene, _) in enumerate(aligned_data):
        genome_name = focal_gene.genome.name if focal_gene.genome else ""
        by_genome[genome_name].append(idx)

    for indices in by_genome.values():
        for a, b in combinations(indices, 2):
            if mode == "tandem":
                should_collapse = _is_tandem(aligned_data[a], aligned_data[b])
            else:
                should_collapse = bool(_window_keys(aligned_data[a][1]) & _window_keys(aligned_data[b][1]))
            if should_collapse:
                union(a, b)

    groups: dict[int, list[int]] = defaultdict(list)
    for idx in range(len(aligned_data)):
        groups[find(idx)].append(idx)

    return [groups[root] for root in sorted(groups)]


def _is_tandem(row_a: tuple[Gene, list[Gene]], row_b: tuple[Gene, list[Gene]]) -> bool:
    """Reports whether two focal genes are drawn next to each other on the same contig.

    Adjacency is judged on the genes actually plotted, not on the genome: genes that sit
    between the two copies but are left out of the figure (for instance genes without an
    orthogroup when `-k/--keep_all_genes` is off) do not break the tandem run.

    Args:
        row_a (tuple[Gene, list[Gene]]): The focal gene and aligned window of the first track.
        row_b (tuple[Gene, list[Gene]]): The focal gene and aligned window of the second track.

    Returns:
        bool: True if no plotted gene of either window lies between the two focal genes.
    """
    gene_a, gene_b = row_a[0], row_b[0]
    if gene_a.seqid != gene_b.seqid or gene_a.index is None or gene_b.index is None:
        return False

    first, last = sorted((gene_a.index, gene_b.index))
    if last - first == 1:
        return True

    plotted = {gene.index for gene in chain(row_a[1], row_b[1]) if gene and gene.index is not None and gene.seqid == gene_a.seqid}
    return not any(first < index < last for index in plotted)


def _window_keys(window: list[Gene]) -> set[tuple[str, int]]:
    """Collects the identity of the genes present in a window.

    Args:
        window (list[Gene]): The genes of an aligned window, padding included.

    Returns:
        set[tuple[str, int]]: One (contig, gene index) key per gene of the window.
    """
    return {_gene_key(gene) for gene in window if gene and gene.index is not None}


def _get_palette(tracks: list[Track]) -> dict[str, str]:
    """Prepares color palette for a SOG.

    Args:
        tracks (list[Track]): The tracks of the figure.

    Returns:
        dict[str, str]: A color palette for this SOG.
    """
    # Count OG occurrences to determine coloring
    og_counter = Counter()
    for track in tracks:
        for gene in track.genes:
            if gene and gene.og:
                og_counter[gene.og.id] += 1

    # Build palette
    ColorCycler.reset_cycle()
    palette = {og_id: ColorCycler() for og_id, count in og_counter.items() if count > 1}

    return palette


def _render_sog_figure(tracks: list[Track], palette: dict[str, str], output_path: Path) -> None:
    """Renders the GenomeViz object and saves the file.

    Args:
        tracks (list[Track]): The tracks of the figure, already shaped by `_build_tracks`.
        palette (dict[str, str]): A color palette for the SOG.
        output_path (Path): The path where the rendered figure should be saved.

    Returns:
        None
    """
    # Get length from the longest track
    length = max(len(track.genes) for track in tracks)

    gv = GenomeViz(fig_width=length * 1.8, track_align_type="left")

    for track in tracks:
        if track.is_separator:
            gv.add_feature_track(track.title, 1)
            continue

        feature_track = gv.add_feature_track(track.title, length * 5 - 2)
        _mark_contig_ends(feature_track, track)

        for idx, gene in enumerate(track.genes):
            if not gene:
                continue

            og_label = gene.og.id if gene.og else "None"
            if og_label == "None":
                fc = "#dddddd"  # Default for genes without OGs
            else:
                fc = palette.get(og_label, "#888888")  # OGs that are not contribute to synteny

            # All focal genes of a collapsed track share the focal gene color
            is_focal = gene.id in track.focal_ids or gene.representative.id in track.focal_ids
            text_weight = "bold" if is_focal else "normal"
            if is_focal:
                fc = _FOCAL_COLOR

            start = idx * 5
            feature_track.add_feature(
                start,
                start + 3,
                plotstyle="bigrbox",
                label=_wrap_label(gene.id),
                fc=fc,
                text_kws={"weight": text_weight, "rotation": 0, "hpos": "center"},
            )
            feature_track.add_text(start + 1.5, og_label, vpos="bottom", hpos="center", rotation=0, weight=text_weight)

    gv.savefig(output_path)


def _mark_contig_ends(feature_track, track: Track) -> None:
    """Draws a bar where the window of a track is cut short by the end of its contig.

    Args:
        feature_track (FeatureTrack): The track being drawn.
        track (Track): The track to read the contig edges from.

    Returns:
        None
    """
    columns = [idx for idx, gene in enumerate(track.genes) if gene]
    if not columns:
        return

    # The margin column reserved by _build_tracks holds the bar
    if track.at_contig_start:
        start = columns[0] * 5 - 2
        feature_track.add_feature(start, start + 1, plotstyle="bigbox", fc=_CONTIG_END_COLOR)

    if track.at_contig_end:
        start = columns[-1] * 5 + 4
        feature_track.add_feature(start, start + 1, plotstyle="bigbox", fc=_CONTIG_END_COLOR)


def _wrap_label(text: str, width: int = _LABEL_WIDTH, min_tail: int = _LABEL_MIN_TAIL) -> str:
    """Wraps a gene label without leaving a stub on the last line.

    `textwrap.wrap` fills every line before moving on, so a label slightly longer than a
    multiple of `width` ends in one or two leftover characters. When that happens the label is
    re-wrapped at the widest smaller width that keeps the same number of lines and gives a
    longer last line, which shifts the break instead of adding a line. Labels that already end
    in a long enough line, such as those broken at a hyphen, are left untouched.

    Args:
        text (str): The label to wrap.
        width (int): The maximum number of characters per line.
        min_tail (int): The shortest acceptable last line.

    Returns:
        str: The label, with its lines separated by newlines.
    """
    lines = textwrap.wrap(text, width=width)
    if len(lines) < 2 or len(lines[-1]) >= min_tail:
        return "\n".join(lines) if lines else text

    for candidate_width in range(width - 1, 0, -1):
        candidate = textwrap.wrap(text, width=candidate_width)
        if len(candidate) == len(lines) and len(candidate[-1]) >= min_tail:
            return "\n".join(candidate)

    # No narrower width helps (e.g. the text has one long unbreakable run)
    return "\n".join(lines)


def _track_label(members: list[Gene], max_shown: int = 3) -> str:
    """Builds the gene part of a track title.

    Args:
        members (list[Gene]): The focal genes the track represents.
        max_shown (int): How many gene IDs to spell out before summarizing the rest.

    Returns:
        str: The gene IDs of the track, one per line, truncated when the track is collapsed
        into more than `max_shown` focal genes.
    """
    gene_ids = [gene.id for gene in members]
    if len(gene_ids) > max_shown:
        gene_ids = gene_ids[:max_shown] + [f"(+{len(members) - max_shown} more)"]
    return "\n".join(gene_ids)


if __name__ == "__main__":
    run_cli()
