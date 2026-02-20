"""
Library module for orthoSynAssign containing file parsing utilities.
"""

from __future__ import annotations

from .gene import Gene, Genome
from .orthogroup import Orthogroup, align_sog_dict, compare_gene_sets, get_shared_ogs
from .parsers import BedParser, read_og_table, write_og_table

__all__ = [
    Gene,
    Genome,
    Orthogroup,
    compare_gene_sets,
    get_shared_ogs,
    align_sog_dict,
    BedParser,
    read_og_table,
    write_og_table,
]
