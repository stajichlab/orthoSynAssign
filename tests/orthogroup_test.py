import pickle

import pytest

from orthosynassign.lib import align_sog_dict, get_shared_ogs


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
        """Test that genes with .og = None are ignored safely."""
        genome_a = genome_factory("Genome_A")
        genome_b = genome_factory("Genome_B")

        og_shared = og_factory("OG1")

        # Gene with an OG
        g_a1 = gene_factory("A1")
        g_b1 = gene_factory("B1")

        # Genes WITHOUT an OG
        g_orphan_a = gene_factory("OrphanA")  # .og is None by default
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


class TestAlignSogDict:
    def test_basic_alignment(self, gene_factory):
        """
        Test that two neighborhoods with different focal gene offsets
        are shifted to match a common pivot.
        """
        g1, g2, g3 = gene_factory("G1"), gene_factory("G2"), gene_factory("G3")
        focal_a, focal_b = gene_factory("focal_a"), gene_factory("focal_B")

        # Scenario:
        # Dict 1: [G1, G2, focal_A] -> focal is at index 2
        # Dict 2: [focal_B, G3]     -> focal is at index 0
        sog_dict = {focal_a: [g1, g2, focal_a], focal_b: [focal_b, g3]}

        aligned = align_sog_dict(sog_dict)

        # The pivot should be 2 (the max index of a focal gene)
        # List 1: [G1, G2, focal_A] (no change to front, needs 0 backpad)
        # List 2: [None, None, focal_B, G3] (needs 2 frontpads to move focal_B to index 2)

        assert aligned[focal_a] == [g1, g2, focal_a, None]
        assert aligned[focal_b] == [None, None, focal_b, g3]
        assert len(aligned[focal_a]) == len(aligned[focal_b])

    def test_empty_dict(self):
        """Ensure it handles an empty input dictionary without crashing."""
        assert align_sog_dict({}) == {}

    def test_focal_gene_not_in_list(self, gene_factory):
        """
        Test fallback logic: if a focal gene isn't in its list,
        it treats its index as 0.
        """
        g_focal = gene_factory("Focal")
        g_other = gene_factory("Other")

        # Focal gene is NOT in the neighborhood list
        sog_dict = {g_focal: [g_other]}

        aligned = align_sog_dict(sog_dict)

        assert aligned[g_focal] == [g_other]

    def test_pure_backpadding(self, gene_factory):
        """
        Tests that lists are padded at the end (back_pad) to
        ensure equal total length.
        """
        focal_a, focal_b = gene_factory("focal_a"), gene_factory("focal_b")
        g1, g2, g3 = gene_factory("G1"), gene_factory("G2"), gene_factory("G3")

        # focal_a is at index 0, length 2: [focal_a, G1]
        # focal_b is at index 0, length 4: [focal_b, G1, G2, G3]
        sog_dict = {focal_a: [focal_a, g1], focal_b: [focal_b, g1, g2, g3]}

        aligned = align_sog_dict(sog_dict)

        # Both should have total length 4
        assert aligned[focal_a] == [focal_a, g1, None, None]
        assert aligned[focal_b] == [focal_b, g1, g2, g3]

    def test_max_offset_at_start(self, gene_factory):
        """Tests alignment when one list is heavily biased to the right."""
        g_far_right = gene_factory("Right")
        g_far_left = gene_factory("Left")
        n = gene_factory("N")

        sog_dict = {
            g_far_right: [n, n, n, g_far_right],  # Index 3
            g_far_left: [g_far_left],  # Index 0
        }

        aligned = align_sog_dict(sog_dict)

        # Left gene list should get 3 Nones in front
        assert aligned[g_far_left] == [None, None, None, g_far_left]
