from unittest.mock import MagicMock

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
    with open("tests/data/orthogroups.tsv") as f:
        header = f.readline().strip().split("\t")
        samples = header[1:]
        og_rows = {}
        for line in f:
            parts = line.strip().split("\t")
            og_id = parts[0]
            og_rows[og_id] = parts[1:]

    genomes = {}
    for sample in samples:
        genome = genome_factory(sample)
        with open(f"tests/data/{sample}.bed") as f:
            for line in f:
                parts = line.strip().split("\t")
                seqid, start, end, gene_id = parts[0], parts[1], parts[2], parts[3]
                genome.add_gene(gene_factory(gene_id, seqid, start, end))
        genomes[sample] = genome

    ogs = []
    for og_id, gene_cols in og_rows.items():
        og = og_factory(og_id)
        for genes, genome in zip(gene_cols, genomes.values()):
            for gene in genes.split(","):
                og.add_gene(genome[gene.strip()])
        ogs.append(og)

    return genomes, ogs
