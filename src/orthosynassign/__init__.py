"""
orthoSynAssign - Ortholog Synteny Assignment Tool

A tool for analyzing orthologous groups and synteny using OrthoFinder or other tool results
and genome annotation files (GFF3/GTF).
"""

from importlib.metadata import metadata

VERSION = metadata("orthosynassign")["Version"]
AUTHOR = metadata("orthosynassign")["Author-email"]
