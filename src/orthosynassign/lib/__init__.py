"""
Library module for orthoSynAssign containing file parsing utilities.
"""

from __future__ import annotations

from .gene import Gene, Genome
from .orthogroup import Orthogroup, align_sog_dict, get_shared_ogs
from .parsers import BedParser, read_og_table, write_og_table
from .synteny import SyntenyEngine, calculate_synteny_ratio, cluster_refined_ogs

__all__ = [
    Gene,
    Genome,
    Orthogroup,
    align_sog_dict,
    get_shared_ogs,
    BedParser,
    read_og_table,
    write_og_table,
    SyntenyEngine,
    calculate_synteny_ratio,
    cluster_refined_ogs,
]
