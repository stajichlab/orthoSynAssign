"""
Library module for orthoSynAssign containing file parsing utilities.
"""

from .parsers import read_gtf, read_gff3, read_orthofinder_table

__all__ = ["read_gtf", "read_gff3", "read_orthofinder_table"]
