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
        assert og[0] == gene
        assert gene.orthogroup == og
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


class TestOrthogroupRefinement:
    def test_get_refined_sogs_integration(self, gene_factory, genome_factory, og_factory, og) -> None:
        """
        Tests the full flow of get_refined_sogs using real functions.
        This ensures Orthogroup, compare_gene_sets, and consolidate_into_sogs
        all talk to each other correctly.
        """
        # 1. Setup: Create two genomes with one perfectly syntenic pair
        genome_a = genome_factory("Genome_A")
        genome_b = genome_factory("Genome_B")

        # We need at least one neighbor to satisfy window_size=2
        # Anchor genes (the focal ones)
        g_a_focal = gene_factory("A_focal", "chr1", 1000, 2000)
        g_b_focal = gene_factory("B_focal", "chr1", 1000, 2000)

        # Syntenic neighbors (to ensure the ratio is 1.0)
        g_a_neighbor = gene_factory("A_neighbor", "chr1", 2100, 3100)
        g_b_neighbor = gene_factory("B_neighbor", "chr1", 2100, 3100)

        # Setup genomic context
        for g in [g_a_focal, g_a_neighbor]:
            genome_a.add_gene(g)
        for g in [g_b_focal, g_b_neighbor]:
            genome_b.add_gene(g)

        # Assign neighbors to a different OG so they act as anchors
        neighbor_og = og_factory("OG_NEIGHBOR")
        neighbor_og.add_gene(g_a_neighbor)
        neighbor_og.add_gene(g_b_neighbor)

        # Add focal genes to the test OG
        og.add_gene(g_a_focal)
        og.add_gene(g_b_focal)

        # 2. Run the actual logic
        # ratio_threshold=1.0 ensures they MUST match perfectly
        result = og.get_refined_sogs(window_size=2, ratio_threshold=1.0)

        # 3. Assertions
        assert len(result) == 1
        sog = result[0]
        assert g_a_focal in sog[genome_a]
        assert g_b_focal in sog[genome_b]

    def test_get_refined_sogs_no_synteny_found(self, gene_factory, genome_factory, og):
        """Test that an OG with no syntenic support returns an empty list."""
        genome_a = genome_factory("Genome_A")
        genome_b = genome_factory("Genome_B")

        # Genes in different scaffolds/locations with no neighbors
        g_a = gene_factory("A1", "chr1", 1000, 2000)
        g_b = gene_factory("B1", "chr2", 5000, 6000)

        genome_a.add_gene(g_a)
        genome_b.add_gene(g_b)
        og.add_gene(g_a)
        og.add_gene(g_b)

        # Since there are no shared neighbors, this should return []
        result = og.get_refined_sogs(window_size=4, ratio_threshold=0.5)
        assert result == []


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
        assert restored[0].orthogroup is restored
