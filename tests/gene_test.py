import pickle

import pytest

# --- Gene ---


class TestGeneInitiation:
    def test_gene_initialization(self, gene_factory):
        """Test that a Gene object initializes correctly with forward coordinates."""
        gene = gene_factory("A1")
        assert gene.seqid == "scaf1"
        assert gene.start == 100
        assert gene.len == 100
        assert gene.id == "A1"
        assert gene.representative == gene

        # Assert pointers default to None
        assert gene.genome is None
        assert gene.og is None
        assert gene.index is None

    def test_gene_assign_representative(self, gene_factory):
        gene = gene_factory("A1")
        isoform = gene_factory("A1_isoform_1")
        isoform.representative = gene
        assert isoform.representative == gene

    def test_gene_reverse_strand_length(self, gene_factory):
        """Test that gene length is strictly positive even if start > end."""
        # Simulating a gene on the reverse/minus strand
        gene = gene_factory("A2_rev", start=2000, end=1000)

        assert gene.start == 2000
        assert gene.len == 1000  # Evaluates abs(1000 - 2000)

    def test_gene_repr_unassigned(self, gene_factory):
        """Test the __repr__ output when genome and orthogroup are None."""
        gene = gene_factory("A1")
        expected_repr = "A1 @ Unknown genome | Unassigned orthogroup"
        assert repr(gene) == expected_repr

    def test_gene_repr_assigned(self, gene_factory, mock_genome_factory, mock_og_factory):
        """Test the __repr__ output when pointers are fully assigned."""
        gene = gene_factory("A1")
        gene.genome = mock_genome_factory("Sample_A")
        gene.og = mock_og_factory("OG001")

        expected_repr = "A1 @ Sample_A | OG001"
        assert repr(gene) == expected_repr

    def test_gene_slots_memory_protection(self, gene_factory):
        """Test that __slots__ prevents dynamic attribute assignment."""
        gene = gene_factory("A1")
        # Gene should not have a __dict__ due to __slots__
        with pytest.raises(AttributeError):
            gene.some_random_attribute = "This should fail"


class TestGeneSerialization:
    def test_gene_pickling(self, gene_factory):
        """Test that __getstate__ and __setstate__ successfully serialize/deserialize the object."""
        original_gene = gene_factory("A1")
        original_gene.index = 42  # Modify a default attribute to ensure it carries over

        # Serialize (simulate passing to a multiprocessing worker)
        pickled_data = pickle.dumps(original_gene)

        # Deserialize
        restored_gene = pickle.loads(pickled_data)

        # Check attributes
        assert restored_gene.seqid == "scaf1"
        assert restored_gene.start == 100
        assert restored_gene.len == 100
        assert restored_gene.id == "A1"
        assert restored_gene.index == 42

        # Pointers should still be None since we didn't mock them for pickling here
        assert restored_gene.genome is None
        assert restored_gene.og is None

    def test_gene_pickling_with_pointers(self, gene_factory, mock_genome_factory, mock_og_factory):
        """Test that genome and orthogroup pointers are preserved during pickling."""
        original_gene = gene_factory("A1")

        # Assign the mock objects to the pointers
        original_gene.genome = mock_genome_factory("Sample_A")
        original_gene.og = mock_og_factory("OG001")

        # Serialize the gene (this will recursively pickle the attached mock objects too)
        pickled_data = pickle.dumps(original_gene)

        # Deserialize the gene
        restored_gene = pickle.loads(pickled_data)

        # 1. Verify the pointers are no longer None
        assert restored_gene.genome is not None
        assert restored_gene.og is not None

        # 2. Verify the data inside the pointed objects survived
        assert restored_gene.genome.name == "Sample_A"
        assert restored_gene.og.id == "OG001"

        # 3. Verify the __repr__ method can successfully read the restored pointers
        assert repr(restored_gene) == "A1 @ Sample_A | OG001"


# --- Genome ---


@pytest.fixture
def empty_genome(genome_factory):
    return genome_factory("Sample_A", chromosome_type="l")


@pytest.fixture
def circular_genome(genome_factory):
    return genome_factory("Sample_Circular", chromosome_type="c")


@pytest.fixture
def populated_genome(gene_factory, empty_genome, mock_og_factory):
    """Creates a linear genome with 10 genes, alternating OG assignments."""
    for i in range(10):
        gene = gene_factory(f"gene_{i}", "chr1", start=i * 100, end=i * 100 + 50)
        if i % 2 == 0:
            gene.og = mock_og_factory(f"OG_{i}")
        empty_genome.add_gene(gene)
    return empty_genome


class TestGenomeInitialization:
    def test_init_valid(self, empty_genome):
        assert empty_genome.name == "Sample_A"
        assert empty_genome.chromosome_type == "l"
        assert len(empty_genome) == 0

    def test_init_invalid_type(self, genome_factory):
        with pytest.raises(ValueError, match="Invalid chromosome type"):
            genome_factory("Fail", chromosome_type="x")


class TestGenomeDataManagement:
    def test_add_gene_and_index_increment(self, gene_factory, empty_genome):
        """Test that multiple genes get unique, incrementing indices."""
        genes = [gene_factory(f"G{i + 1}", start=100 + (i * 200), end=200 + (i * 200)) for i in range(3)]

        for i, gene in enumerate(genes):
            empty_genome.add_gene(gene)
            assert gene.index == i  # Check current index
            assert len(empty_genome) == i + 1

        assert empty_genome[2].id == "G3"

    def test_add_gene_with_isoform(self, gene_factory, empty_genome):
        """Test adding genes with is_isoform tags."""
        genes = [gene_factory(f"G{i + 1}", start=100 + (i * 200), end=200 + (i * 200)) for i in range(3)]

        for i, gene in enumerate(genes):
            empty_genome.add_gene(gene)
        genome_length = len(empty_genome)

        for i in range(3):
            isoform = gene_factory(f"G2_isoform_{i + 1}", start=300, end=400)
            isoform.representative = genes[1]
            empty_genome.add_gene(isoform, is_isoform=True)
            assert isoform.index is None
            assert len(empty_genome) == genome_length
            assert isoform.representative.index == 1

    def test_getitem_retrieval_types(self, populated_genome):
        """Test retrieving genes via index, string ID, and slices."""
        # 1. Retrieve by index (int)
        gene_by_idx = populated_genome[0]
        assert gene_by_idx.id == "gene_0"

        # 2. Retrieve by ID (str)
        gene_by_id = populated_genome["gene_5"]
        assert gene_by_id.id == "gene_5"
        assert gene_by_id.index == 5

        # 3. Retrieve by slice
        # Get first 3 genes
        subset = populated_genome[0:3]
        assert isinstance(subset, list)
        assert len(subset) == 3
        assert [g.id for g in subset] == ["gene_0", "gene_1", "gene_2"]

        # Get every second gene
        stride_subset = populated_genome[::2]
        assert len(stride_subset) == 5
        assert stride_subset[1].id == "gene_2"

    def test_getitem_keyerror(self, empty_genome):
        """Ensure KeyError is raised for missing string IDs."""
        with pytest.raises(KeyError):
            _ = empty_genome["NonExistent"]

    def test_getitem_indexerror(self, populated_genome):
        """Ensure IndexError is raised for out-of-bounds integer access."""
        with pytest.raises(IndexError):
            _ = populated_genome[99]

    def test_iteration(self, populated_genome):
        """Verify the genome is iterable and returns genes in order."""
        for i, gene in enumerate(populated_genome):
            assert gene.id == f"gene_{i}"
            assert gene.index == i


class TestGenomeGetWindow:
    def test_get_window_linear_boundaries(self, populated_genome):
        """Test that linear genomes stop at the start/end of the gene list."""
        focal = populated_genome[0]  # gene_0
        target_ogs = {f"OG_{i}" for i in range(10)}

        # Request window of 4, but at index 0, there is only 'downstream'
        window = populated_genome.get_window(focal, target_ogs, window_size=4)

        # Should find OG_2 and OG_4 (the next two anchored genes)
        assert len(window) == 2
        assert [g.id for g in window] == ["gene_2", "gene_4"]

    def test_get_window_scaffold_break(self, gene_factory, empty_genome, mock_og_factory):
        """Test that windows stop when crossing a scaffold boundary."""
        g1 = gene_factory("G1", "scaf_1", 100, 200)
        g1.og = mock_og_factory("OG1")

        g2 = gene_factory("G2", "scaf_2", 100, 200)  # Different scaffold
        g2.og = mock_og_factory("OG2")

        empty_genome.add_gene(g1)
        empty_genome.add_gene(g2)

        # We look for OG2, but it's on scaf_2, so it should be blocked
        window = empty_genome.get_window(g1, {"OG2"}, window_size=2)
        assert len(window) == 0

    def test_get_window_circularity(self, gene_factory, circular_genome, mock_og_factory):
        """Test that circular genomes wrap around using modulo."""
        # Create 5 genes: G0, G1, G2, G3, G4
        for i in range(5):
            gene = gene_factory(f"G{i}", "chr1", i * 100, i * 100 + 50)
            gene.og = mock_og_factory(f"OG{i}")
            circular_genome.add_gene(gene)

        focal = circular_genome[0]  # G0
        target_ogs = {"OG1", "OG4"}  # OG1 is index 1, OG4 is index 4 (upstream wrap)

        window = circular_genome.get_window(focal, target_ogs, window_size=2)

        # In a circular genome:
        # Downstream 1 step: G1 (OG1 found)
        # Upstream 1 step: G4 (OG4 found)
        ids = [g.id for g in window]
        assert "G4" in ids
        assert "G1" in ids
        assert len(window) == 2

    def test_get_window_isoforms(self, populated_genome, gene_factory):
        """Test get_window using isoform."""
        focal = populated_genome[0]  # gene_0
        isoforms = [f"isoform_{i + 1}" for i in range(3)]
        target_ogs = {f"OG_{i}" for i in range(10)}
        for isoform in isoforms:
            g = gene_factory(isoform)
            g.representative = focal
            populated_genome.add_gene(g, is_isoform=True)

            # Request window of 4
            window = populated_genome.get_window(g, target_ogs, window_size=4)

            # Should find OG_2 and OG_4 (the next two anchored genes)
            assert len(window) == 2
            assert [g.id for g in window] == ["gene_2", "gene_4"]


class TestGenomeSerialization:
    def test_genome_pickle(self, populated_genome):
        """Ensure back-pointers are restored after pickling."""
        data = pickle.dumps(populated_genome)
        restored = pickle.loads(data)

        assert restored.name == populated_genome.name
        assert len(restored) == 10
        # Check back-pointer logic in __setstate__
        assert restored[0].genome == restored
