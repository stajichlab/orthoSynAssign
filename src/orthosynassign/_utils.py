"""
Orthosynassign internal utilities and decorators for use within the package.
"""

from __future__ import annotations

import argparse
import logging
import re
from pathlib import Path
from typing import TYPE_CHECKING, Protocol

from .lib import BedParser

if TYPE_CHECKING:
    from .lib.parsers import AnnotationParser

logger = logging.getLogger(__name__)


class _BaseArgs(Protocol):
    """Arguments shared by ALL entry points.

    Attributes:
        og_file (Path): The Orthogroup file path.
        bed (list[Path]): A list of BED files paths.
        output (Path): The output file path.
        window (int): The size of the genomic window to consider around each orthologous group.
        verbose (bool): Flag to enable verbose logging.
    """

    og_file: Path
    bed: list[Path]
    output: Path
    window: int
    verbose: bool


class RefineArgs(_BaseArgs, Protocol):
    """Arguments specific to refine.

    Attributes:
        threshold (float): The confidence threshold for filtering orthologous groups.
        threads (int): Number of threads to use for parallel processing.
    """

    threshold: float
    threads: int


class VisualizeArgs(_BaseArgs, Protocol):
    """Arguments specific to visualize.

    Attributes:
        sog_file (Path): The Sog file path.
        sog (str): The Sog identifier.
        fmt (str): The format for the visualization output.
        keep_all_genes (bool): Flag to keep all genes in the visualization output.
    """

    sog_file: Path
    sog: str
    fmt: str
    keep_all_genes: bool


class CustomHelpFormatter(argparse.RawDescriptionHelpFormatter):
    """Custom help formatter for argparse to enhance text wrapping and default value display.

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


def setup_logging(verbose: bool = False) -> None:
    """Configure logging for the application.

    Args:
        verbose (bool): Flag to enable verbose logging.
    """
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def validate_annotations(args: _BaseArgs) -> list[AnnotationParser]:
    """
    Validate input files and directories.

    Args:
        args (_BaseArgs): Parsed command line arguments

    Returns:
        list[AnnotationParser]: A list of AnnotationParser objects.
    """
    # Check annotation files
    if getattr(args, "bed", None):
        files = args.bed
        parser = BedParser

    annotations = []
    for file in files:
        logger.debug("Reading annotation from: %s", file)
        annotations.append(parser(file))

    return annotations


def validate_orthogroup(og_input: Path) -> Path:
    """Validate the OrthoFinder file.

    Args:
        og_input (Path): The Orthogroup file path to be validated.

    Returns:
        Path: The validated Orthogroup file path.

    Raises:
        FileNotFoundError: If the Orthogroup file does not exist.
        ValueError: If the Orthogroup file is not a regular file.
    """
    # Check OrthoFinder file
    file = Path(og_input)
    if not file.exists():
        raise FileNotFoundError(f"Orthogroup file not found: {file}")
    if not file.is_file():
        raise ValueError(f"Orthogroup file must be a regular file: {file}")
    return file
