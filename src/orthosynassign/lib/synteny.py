from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from .rs import SyntenyEngine

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from .gene import Genome
    from .orthogroup import Orthogroup


def prepare_and_init_engine(genomes: list[Genome], orthogroups: list[Orthogroup]) -> SyntenyEngine:
    """
    Converts biological objects into integer vectors and initializes
    the Rust SyntenyEngine.
    """
    # 1. Map OG IDs to integers for the entire project
    og_str_to_int = {og.id: i for i, og in enumerate(orthogroups)}

    og_arr_list: list[NDArray[np.int32]] = []
    seq_arr_list: list[NDArray[np.int16]] = []

    for genome in genomes:
        # Local map for scaffold IDs per genome (e.g., 'Chr1' -> 0)
        seqid_map = {}
        next_seqid_int = 0

        og_arr = []
        seq_arr = []

        for gene in genome._genes:
            # Handle scaffold/sequence mapping
            seqid = gene.seqid
            if seqid not in seqid_map:
                seqid_map[seqid] = next_seqid_int
                next_seqid_int += 1
            seq_arr.append(seqid_map[seqid])

            # Handle Orthogroup mapping
            # Use -1 for genes not assigned to an orthogroup
            og_id = gene.og.id if gene.og else None
            og_arr.append(og_str_to_int.get(og_id, -1))

        # We pass lists to the Rust constructor; Rust converts them to Vec internally.
        # Ensure types match: og_list (int32), seq_list (int16)
        og_arr_list.append(np.array(og_arr, dtype=np.int32))
        seq_arr_list.append(np.array(seq_arr, dtype=np.int16))

    # Initialize the Rust Engine
    engine = SyntenyEngine(len(genomes), len(orthogroups), og_arr_list, seq_arr_list)

    return engine
