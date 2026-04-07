import pickle

import pytest


@pytest.fixture
def og(og_factory):
    """Provides a fresh Orthogroup instance."""
    return og_factory("OG00001")


# --- Orthogroup ---


class TestOrthogroupBasics:
    def test_initialization(self, og):
        assert og.id == "OG00001"
        assert len(og) == 0

    def test_add_gene(self, og, gene_factory, genome_factory):
        # 1. Real objects setup
        genome = genome_factory("Genome_A")
        gene = gene_factory("G1")

        # 2. Add to genome first (sets index/genome pointer)
        genome.add_gene(gene)

        # 3. Add to orthogroup
        og.add_gene(gene)

        assert len(og) == 1
        assert gene in og
        assert og[0] == gene
        assert gene.og == og
        assert gene.genome == genome
        assert gene.index == 0

    def test_add_duplicate_gene(self, og, gene_factory):
        gene = gene_factory("G1")
        og.add_gene(gene)
        og.add_gene(gene)  # Should be ignored by the 'if gene_obj not in self._genes' logic
        assert len(og) == 1

    def test_iteration(self, og, gene_factory, genome_factory):
        genome = genome_factory("Genome_A")
        genes = []
        for i in range(3):
            g = gene_factory(f"G{i}")
            genome.add_gene(g)
            og.add_gene(g)
            genes.append(g)

        assert list(og) == genes


class TestOrthogroupSerialization:
    def test_pickle_consistency(self, og, gene_factory, genome_factory):
        """Verify that Orthogroup data and gene associations survive pickling."""
        genome = genome_factory("Genome_A")
        g1 = gene_factory("G1")
        genome.add_gene(g1)
        og.add_gene(g1)

        pickled = pickle.dumps(og)
        restored = pickle.loads(pickled)

        assert restored.id == og.id
        assert len(restored) == 1
        assert restored[0].id == "G1"
        # Ensure the back-pointer to the Orthogroup is still the restored object
        assert restored[0].og is restored


@pytest.fixture
def prepared_comparison(request, read_example_files):
    """
    Bridge fixture: Converts lists of string IDs into lists of Gene objects.
    """
    genome_names, gene_id_lists, expected_id_list = request.param
    genomes, _ = read_example_files

    # Resolve Genome objects
    genome_a = genomes[genome_names[0]]
    genome_b = genomes[genome_names[1]]

    # Resolve Focal Gene Lists (handling potential missing genes)
    genes_a = [genome_a[gid] for gid in gene_id_lists[0]]
    genes_b = [genome_b[gid] for gid in gene_id_lists[1]]

    # Resolve Expected Result List
    # Returns a list of tuples (GeneA, GeneB) if expected exists, else empty list
    expected = []
    if expected_id_list:
        # Assuming the expected_id_list is formatted as ([idA], [idB])
        for id_a, id_b in zip(expected_id_list[0], expected_id_list[1]):
            expected.append((genome_a[id_a], genome_b[id_b]))

    return (genes_a, genes_b), expected
