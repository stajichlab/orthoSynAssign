"""
Orthology utilities for working with GFF features.
"""

from typing import Dict, List, Tuple

from .parsers import GFFFeature


def organize_chromosomes(
    features: List[GFFFeature],
) -> Tuple[Dict[str, int], List[GFFFeature]]:
    """
    Filter gene features, sort by chromosome and start, and index by gene name.

    Args:
        features: Collection of GFFFeature objects to process.

    Returns:
        A tuple (genomicindex, chromsorted) where:
        - genomicindex maps gene name to its index in the sorted array.
        - chromsorted is the list of gene features sorted by seqid then start.
    """
    # this may not work for GTF files we will need to see if that is required or if we will impute those
    gene_features = [feat for feat in features if feat.type == "gene"]

    # Sort genes by chromosome (seqid) then genomic start coordinate
    chromsorted = sorted(gene_features, key=lambda feat: (feat.seqid, feat.start))

    genomicindex: Dict[str, int] = {}
    for idx, feat in enumerate(chromsorted):
        gene_name = (
            feat.attributes.get("Name")
            or feat.attributes.get("gene_id")
            or feat.attributes.get("ID")
        )
        if gene_name:
            genomicindex[gene_name] = idx

    return genomicindex, chromsorted
