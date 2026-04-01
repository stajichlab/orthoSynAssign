import numpy as np
import pytest

from orthosynassign.lib import calculate_synteny_ratio, prepare_and_init_engine


class TestCalculateSyntenyRatio:
    @pytest.mark.parametrize(
        "win_a, win_b, expected_ratio, description",
        [
            # 1. Perfect Identity
            ([1, 2], [1, 2], 1.0, "Identical windows"),
            # 2. Partial Overlap
            ([1, 2, 3], [1, 2, 4], 2 / 3, "Partial overlap (2/3)"),
            # 3. Tandem Duplication (1-to-1 matching)
            ([1, 1], [1], 0.5, "A has extra; match is 1/2"),
            ([1, 1], [1, 1], 1.0, "Both have two copies; match is 2/2"),
            # 4. Order Independence
            ([1, 2], [2, 1], 1.0, "Different order, same content"),
            # 5. Length Penalization
            ([1], [1, 2, 3, 4], 0.25, "B is much longer; ratio drops"),
            # 6. No Overlap
            ([1, 2], [3, 4], 0.0, "Zero shared orthogroups"),
            # 7. Empty Inputs
            ([], [1], 0.0, "Window A is empty"),
            ([1], [], 0.0, "Window B is empty"),
            ([], [], 0.0, "Both windows are empty"),
        ],
    )
    def test_calculate_synteny_ratio(self, win_a, win_b, expected_ratio, description):
        """
        Tests the synteny ratio calculation across multiple genomic scenarios.
        Using pytest.approx for floating point comparisons.
        """
        arr_a = np.array(win_a, dtype=np.int32)
        arr_b = np.array(win_b, dtype=np.int32)

        result = calculate_synteny_ratio(arr_a, arr_b)
        assert result == pytest.approx(expected_ratio), f"Failed: {description}"


class TestSyntenyEngineRefinement:
    @pytest.fixture
    def og(self, og_factory):
        """Provides a fresh Orthogroup instance."""
        return og_factory("OG00001")

    def test_refine_integration(self, gene_factory, genome_factory, og_factory, og) -> None:
        """
        Tests the full flow of refine using real functions.
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

        # Add focal genes to the test OG
        og.add_gene(g_a_focal)
        og.add_gene(g_b_focal)

        # Assign neighbors to a different OG so they act as anchors
        neighbor_og = og_factory("OG_NEIGHBOR")
        neighbor_og.add_gene(g_a_neighbor)
        neighbor_og.add_gene(g_b_neighbor)

        # Initialize the Engine
        # The engine needs the list of all relevant genomes and OGs
        genomes = [genome_a, genome_b]
        orthogroups = [og, neighbor_og]
        engine = prepare_and_init_engine(genomes, orthogroups)
        result = engine.refine(0, window_size=2, ratio_threshold=1.0)

        # 3. Assertions
        assert len(result) == 1
        sog = result[0]
        sog_genes = [genomes[genome_idx][gene_idx] for genome_idx, gene_idx in sog]
        assert g_a_focal in sog_genes
        assert g_b_focal in sog_genes

    def test_refine_no_synteny_found(self, gene_factory, genome_factory, og):
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
        genomes = [genome_a, genome_b]
        orthogroups = [og]
        engine = prepare_and_init_engine(genomes, orthogroups)
        result = engine.refine(0, window_size=2, ratio_threshold=1.0)
        assert result == []
