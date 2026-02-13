import argparse
import logging
import re
from pathlib import Path

from .lib.parsers import AnnotationParser, BedParser


def setup_logging(verbose: bool = False):
    """Configure logging for the application."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


class CustomHelpFormatter(argparse.RawDescriptionHelpFormatter):
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


def validate_annotations(args, logger) -> list[AnnotationParser]:
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


def validate_orthogroup(args) -> Path:
    # Check OrthoFinder file
    file = Path(args.og_input)
    if not file.exists():
        raise FileNotFoundError(f"Orthogroup file not found: {file}")
    if not file.is_file():
        raise ValueError(f"Orthogroup file must be a regular file: {file}")
    return file
