import pytest

from orthosynassign.lib import Gene, Genome, Orthogroup


class MockGenome:
    def __init__(self, name: str):
        self.name = name


class MockOrthogroup:
    def __init__(self, og_id: str):
        self.id = og_id


@pytest.fixture
def mock_genome_factory():
    """Returns a function that creates MockGenome objects on demand."""

    def _make(name: str):
        return MockGenome(name)

    return _make


@pytest.fixture
def mock_og_factory():
    """Returns a function that creates MockOrthogroup objects on demand."""

    def _make(og_id: str):
        return MockOrthogroup(og_id)

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

    def _make(name: str, chromosome_type="l"):
        return Genome(name, chromosome_type=chromosome_type)

    return _make


@pytest.fixture
def og_factory():
    """Creates a real Orthogroup instance."""

    def _make(og_id: str):
        return Orthogroup(og_id)

    return _make
