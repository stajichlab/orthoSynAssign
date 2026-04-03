from __future__ import annotations

from typing import TYPE_CHECKING

from .rs import SyntenyEngine

if TYPE_CHECKING:
    from .gene import Genome
    from .orthogroup import Orthogroup


def prepare_and_init_engine(genomes: list[Genome], orthogroups: list[Orthogroup]) -> SyntenyEngine:
    """
    Converts biological objects into integer vectors and initializes
    the Rust SyntenyEngine.
    """
    # 1. Map OG IDs to integers for the entire project
    og_str_to_int = {og.id: i for i, og in enumerate(orthogroups)}

    og_list_all: list[list[int]] = []
    seqid_list_all: list[list[int]] = []

    for genome in genomes:
        # Local map for scaffold IDs per genome (e.g., 'Chr1' -> 0)
        seqid_map: dict[str, int] = {}
        next_seqid_int = 0

        og_list_genome: list[int] = []
        seqid_list_genome: list[int] = []

        for gene in genome._genes:
            # Handle scaffold/sequence mapping
            seqid = gene.seqid
            if seqid not in seqid_map:
                seqid_map[seqid] = next_seqid_int
                next_seqid_int += 1
            seqid_list_genome.append(seqid_map[seqid])

            # Handle Orthogroup mapping
            # Use -1 for genes not assigned to an orthogroup
            og_id = gene.og.id if gene.og else None
            og_list_genome.append(og_str_to_int.get(og_id, -1))

        # We pass lists to the Rust constructor; Rust converts them to Vec internally.
        og_list_all.append(og_list_genome)
        seqid_list_all.append(seqid_list_genome)

    # Initialize the Rust Engine
    engine = SyntenyEngine(len(orthogroups), og_list_all, seqid_list_all)

    return engine
