import pickle

import pytest

from orthosynassign.lib import compare_gene_sets, get_shared_ogs


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


class TestCompareGeneSets:
    @pytest.mark.parametrize(
        "prepared_comparison",
        [
            # Case 1: Perfect match (A6-B6)
            # Format: ((Genomes), ([ListA], [ListB]), ([ExpectedA], [ExpectedB]))
            (("Sample_A", "Sample_B"), (["A6"], ["B6"]), (["A6"], ["B6"])),
            # Case 2: No match (Noise)
            (("Sample_A", "Sample_B"), (["A14"], ["B14"]), None),
            # Case 3: Multiple candidates (Paralog selection)
            # Testing if A2 picks the correct B candidate from a list
            (("Sample_A", "Sample_B"), (["A2"], ["B2", "B2b"]), (["A2"], ["B2"])),
            # Case 4: Multiple candidates (Entire duplication)
            (("Sample_A", "Sample_B"), (["A23", "A23b"], ["B23", "B23b"]), (["A23", "A23b"], ["B23", "B23b"])),
        ],
        indirect=["prepared_comparison"],
    )
    def test_compare_gene_sets_using_example_files(self, prepared_comparison):
        """Test syntenic identification with list-based inputs."""
        (genes_a, genes_b), expected = prepared_comparison

        # window_size and ratio_threshold match your example file parameters
        refined_results = compare_gene_sets(genes_a, genes_b, window_size=4, ratio_threshold=0.5)

        # We compare the list of tuples directly
        assert refined_results == expected


class TestGetSharedOGs:
    def test_basic_intersection(self, genome_factory, gene_factory, og_factory):
        """Test that only OGs present in both genomes are returned."""
        genome_a = genome_factory("Genome_A")
        genome_b = genome_factory("Genome_B")

        # Shared OGs
        og_shared = og_factory("OG_SHARED")
        # Unique OGs
        og_a_only = og_factory("OG_A")
        og_b_only = og_factory("OG_B")

        # Setup Genome A: has Shared and A_only
        g1 = gene_factory("G1")
        g2 = gene_factory("G2")
        genome_a.add_gene(g1)
        genome_a.add_gene(g2)

        # Setup Genome B: has Shared and B_only
        g3 = gene_factory("G3")
        g4 = gene_factory("G4")
        genome_b.add_gene(g3)
        genome_b.add_gene(g4)

        og_shared.add_gene(g1)
        og_a_only.add_gene(g2)
        og_shared.add_gene(g3)
        og_b_only.add_gene(g4)

        shared = get_shared_ogs(genome_a, genome_b)

        assert shared == {"OG_SHARED"}
        assert "OG_A" not in shared
        assert "OG_B" not in shared

    def test_handling_missing_orthogroups(self, genome_factory, gene_factory, og_factory):
        """Test that genes with .orthogroup = None are ignored safely."""
        genome_a = genome_factory("Genome_A")
        genome_b = genome_factory("Genome_B")

        og_shared = og_factory("OG1")

        # Gene with an OG
        g_a1 = gene_factory("A1")
        g_b1 = gene_factory("B1")

        # Genes WITHOUT an OG
        g_orphan_a = gene_factory("OrphanA")  # .orthogroup is None by default
        g_orphan_b = gene_factory("OrphanB")

        genome_a.add_gene(g_a1)
        genome_a.add_gene(g_orphan_a)
        genome_b.add_gene(g_b1)
        genome_b.add_gene(g_orphan_b)

        og_shared.add_gene(g_a1)
        og_shared.add_gene(g_b1)

        shared = get_shared_ogs(genome_a, genome_b)

        # Should only contain OG1, and should not crash/contain None
        assert shared == {"OG1"}

    def test_no_overlap(self, genome_factory, gene_factory, og_factory):
        """Test return value is an empty set when no OGs are shared."""
        genome_a = genome_factory("Genome_A")
        genome_b = genome_factory("Genome_B")
        og_a = og_factory("OG_A")
        og_b = og_factory("OG_B")

        g1 = gene_factory("A")
        g2 = gene_factory("B")
        genome_a.add_gene(g1)
        genome_b.add_gene(g2)

        og_a.add_gene(g1)
        og_b.add_gene(g2)

        assert get_shared_ogs(genome_a, genome_b) == set()

    def test_empty_genomes(self, genome_factory):
        """Test that empty genomes return an empty set."""
        genome_a = genome_factory("Empty_A")
        genome_b = genome_factory("Empty_B")

        assert get_shared_ogs(genome_a, genome_b) == set()
