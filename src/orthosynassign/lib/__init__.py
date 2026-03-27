"""
Library module for orthoSynAssign containing file parsing utilities.
"""

from __future__ import annotations

from .gene import Gene, Genome
from .orthogroup import SOG, Orthogroup, align_sog_dict, calculate_synteny_ratio, compare_gene_sets, get_shared_ogs
from .parsers import BedParser, read_og_table, write_og_table
from .synteny import SyntenyEngine

__all__ = [
    Gene,
    Genome,
    Orthogroup,
    SOG,
    compare_gene_sets,
    get_shared_ogs,
    align_sog_dict,
    calculate_synteny_ratio,
    BedParser,
    read_og_table,
    write_og_table,
    SyntenyEngine,
]
