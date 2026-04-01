"""
Library module for orthoSynAssign containing file parsing utilities.
"""

from __future__ import annotations

from .gene import Gene, Genome
from .orthogroup import Orthogroup, align_sog_dict, get_shared_ogs
from .parsers import BedParser, read_og_table, write_og_table
from .rs import calculate_synteny_ratio
from .synteny import prepare_and_init_engine

__all__ = [
    Gene,
    Genome,
    Orthogroup,
    align_sog_dict,
    get_shared_ogs,
    BedParser,
    read_og_table,
    write_og_table,
    calculate_synteny_ratio,
    prepare_and_init_engine,
]
