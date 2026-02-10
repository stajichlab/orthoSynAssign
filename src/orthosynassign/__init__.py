"""
orthoSynAssign - Ortholog Synteny Assignment Tool

Refine OrthoFinder orthologous groups with synteny analysis using genome annotations.
"""

import logging
from importlib.metadata import metadata

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

VERSION = metadata("orthosynassign")["Version"]
AUTHOR = metadata("orthosynassign")["Author-email"]
