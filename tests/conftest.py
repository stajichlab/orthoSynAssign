from unittest.mock import MagicMock

import pandas as pd
import pytest

from orthosynassign.lib import Gene, Genome, Orthogroup


@pytest.fixture
def mock_genome_factory():
    """Returns a function that creates MockGenome objects on demand."""

    def _make(name: str):
        mock_cls = MagicMock(spec=Genome)
        mock_obj = MagicMock(return_value=mock_cls)
        mock_obj.name = name
        return mock_obj

    return _make


@pytest.fixture
def mock_og_factory():
    """Returns a function that creates MockOrthogroup objects on demand."""

    def _make(og_id: str):
        mock_cls = MagicMock(spec=Orthogroup)
        mock_obj = MagicMock(return_value=mock_cls)
        mock_obj.id = og_id
        return mock_obj

    return _make


@pytest.fixture
def gene_factory():
    """Creates a real Gene instance."""

    def _make(gene_id: str, seqid="scaf1", start=100, end=200):
        return Gene(seqid, start, end, gene_id)

    return _make


@pytest.fixture
def genome_factory():
    """Creates a real Genome instance."""

    def _make(name: str, is_circular=False):
        return Genome(name, is_circular=is_circular)

    return _make


@pytest.fixture
def og_factory():
    """Creates a real Orthogroup instance."""

    def _make(og_id: str):
        return Orthogroup(og_id)

    return _make


@pytest.fixture
def read_example_files(gene_factory, genome_factory, og_factory):
    ogs_df: pd.DataFrame = pd.read_csv("tests/data/orthogroups.tsv", sep="\t", index_col=0)
    genomes = {}
    for sample in ogs_df.columns:
        df: pd.DataFrame = pd.read_csv(f"tests/data/{sample}.bed", sep="\t", header=None)
        df = df[[3, 0, 1, 2]]
        genome = genome_factory(sample)
        for idx, row in df.iterrows():
            genome.add_gene(gene_factory(*row.values))
        genomes[sample] = genome

    ogs = []
    for og_id, row in ogs_df.iterrows():
        og = og_factory(og_id)
        for genes, genome in zip(row.values, genomes.values()):
            for gene in genes.split(","):
                og.add_gene(genome[gene.strip()])
        ogs.append(og)

    return genomes, ogs
