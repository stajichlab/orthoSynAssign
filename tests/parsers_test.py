import gzip

import pytest

from orthosynassign.lib import BedParser, read_og_table, write_og_table


@pytest.fixture
def mock_bed_content():
    """Returns a standard 4-column BED string."""
    return "chr1\t100\t200\tGeneA\n# This is a comment\n\nchr1\t300\t400\tGeneB\nchr2\t500\t600\tGeneC\n"


class TestBedParser:
    def test_init_raises_errors(self, tmp_path):
        """Verify that missing files or directories raise appropriate errors."""
        non_existent = tmp_path / "missing.bed"
        directory = tmp_path / "folder"
        directory.mkdir()

        with pytest.raises(FileNotFoundError):
            BedParser(non_existent)

        with pytest.raises(ValueError, match="must be a regular file"):
            BedParser(directory)

    def test_parse_plain_text(self, tmp_path, mock_bed_content):
        """Verify parsing of a standard uncompressed BED file."""
        bed_file = tmp_path / "sample.bed"
        bed_file.write_text(mock_bed_content)

        parser = BedParser(bed_file)
        genome = parser.parse()

        assert len(genome) == 3
        assert genome.name == "sample"
        assert genome[0].id == "GeneA"
        assert genome[2].seqid == "chr2"

    def test_parse_gzipped(self, tmp_path, mock_bed_content):
        """Verify the gzip detection and transparent decompression."""
        bed_file = tmp_path / "sample.bed.gz"
        with gzip.open(bed_file, "wt") as f:
            f.write(mock_bed_content)

        parser = BedParser(bed_file)
        # Check internal magic number detection
        assert parser._is_gzip_file() is True

        genome = parser.parse()
        assert len(genome) == 3
        assert genome[1].id == "GeneB"

    def test_parser_invalid_format(self, tmp_path):
        """Verify that incorrect number of fields raises ValueError."""
        bad_content = "chr1\t100\t200\tGeneA\tExtraField\n"
        bed_file = tmp_path / "invalid.bed"
        bed_file.write_text(bad_content)

        parser = BedParser(bed_file)
        with pytest.raises(ValueError, match="Expected 4 fields"):
            parser.parse()

    def test_parse_exception_resets_genome(self, tmp_path):
        """Verify that if parsing fails, the internal _genome is reset to None."""
        bed_file = tmp_path / "fail.bed"
        bed_file.write_text("not\ta\tbed\tfile\textra")

        parser = BedParser(bed_file)
        try:
            parser.parse()
        except Exception:
            pass

        assert parser._genome is None

    def test_comments_and_empty_lines(self, tmp_path):
        """Ensure comments (#) and blank lines are skipped correctly."""
        content = "\n# header\nchr1\t10\t20\tG1\n\n   \nchr1\t30\t40\tG2"
        bed_file = tmp_path / "clean.bed"
        bed_file.write_text(content)

        parser = BedParser(bed_file)
        genome = parser.parse()

        assert len(genome) == 2
        assert [g.id for g in genome] == ["G1", "G2"]


class TestReadOGTable:
    @pytest.fixture
    def mock_genomes(self, genome_factory, gene_factory):
        """Creates a mock environment with two genomes and pre-loaded genes."""
        ga = genome_factory("Sample_A")
        gb = genome_factory("Sample_B")

        # Add genes to the internal map so the parser can find them
        ga.add_gene(gene_factory("geneA1", "chr1", 10, 20))
        ga.add_gene(gene_factory("geneA2", "chr1", 30, 40))
        gb.add_gene(gene_factory("geneB1", "chr1", 10, 20))

        return [ga, gb]

    def test_standard_parsing(self, tmp_path, mock_genomes):
        """Verify parsing of a standard TSV with multiple genes per cell."""
        og_file = tmp_path / "orthogroups.tsv"
        # Note the mixed delimiters: comma and space
        og_file.write_text("Orthogroup\tSample_A\tSample_B\nOG001\tgeneA1, geneA2\tgeneB1\n")

        ogs = read_og_table(og_file, mock_genomes)

        assert len(ogs) == 1
        og = ogs[0]
        assert og.id == "OG001"
        assert len(og) == 3
        # Verify gene pointer was set
        assert all(g.og == og for g in og)

    def test_missing_sample_error(self, tmp_path, mock_genomes):
        """Verify that a sample in the TSV not in genomes dict raises ValueError."""
        og_file = tmp_path / "orthogroups.tsv"
        og_file.write_text("Orthogroup\tSample_A\tSample_Unknown\nOG001\tgeneA1\tgeneX\n")

        with pytest.raises(ValueError, match="Samples in OG file not found"):
            read_og_table(og_file, mock_genomes)

    def test_gene_not_found_warning(self, tmp_path, mock_genomes, caplog):
        """Verify that if a gene ID isn't in the Genome object, it logs a warning."""
        og_file = tmp_path / "orthogroups.tsv"
        og_file.write_text("Orthogroup\tSample_A\tSample_B\nOG001\tMissingGene\tgeneB1\n")

        ogs = read_og_table(og_file, mock_genomes)

        assert "Gene MissingGene not found in Sample_A" in caplog.text
        assert len(ogs[0]) == 1  # Only geneB1 should be added

    def test_empty_cells(self, tmp_path, mock_genomes):
        """Verify that empty cells (missing genes in a sample) are handled."""
        og_file = tmp_path / "orthogroups.tsv"
        og_file.write_text("Orthogroup\tSample_A\tSample_B\nOG001\tgeneA1\t\nOG002\t\tgeneB1\n")

        ogs = read_og_table(og_file, mock_genomes)
        assert len(ogs) == 2
        assert len(ogs[0]) == 1  # Only Sample_A gene
        assert len(ogs[1]) == 1  # Only Sample_B gene


class TestWriteOGTable:
    def test_write_og_table_output(self, tmp_path, gene_factory, genome_factory):
        """
        Verify the TSV structure, headers, and paralog comma-separation.
        """
        # 1. Setup Mock Data
        ga = genome_factory("Sample_A")
        gb = genome_factory("Sample_B")
        gc = genome_factory("Sample_C")

        # SOG 1: Simple 1-to-1-to-1
        g1a = gene_factory("G1A")
        ga.add_gene(g1a)
        g1b = gene_factory("G1B")
        gb.add_gene(g1b)
        g1c = gene_factory("G1C")
        gc.add_gene(g1c)
        sog1 = ("SOG001", [g1a, g1b, g1c])

        # SOG 2: Paralog in A, missing in C
        g2a1 = gene_factory("G2A1")
        ga.add_gene(g2a1)
        g2a2 = gene_factory("G2A2")
        ga.add_gene(g2a2)
        g2b = gene_factory("G2B")
        gb.add_gene(g2b)
        sog2 = ("SOG002", [g2a1, g2a2, g2b])

        # 2. Run the function
        output_file = tmp_path / "Refined_SOGs.tsv"
        results_iterator = iter([sog1, sog2])
        genome_names = ["Sample_A", "Sample_B", "Sample_C"]

        write_og_table(results_iterator, genome_names, output_file)

        # 3. Assertions
        assert output_file.exists()
        lines = output_file.read_text().splitlines()

        # Check Header
        assert lines[0] == "SOG_ID\tSample_A\tSample_B\tSample_C"

        # Check Row 1 (1-to-1-to-1)
        # Note: Gene IDs are used in the join
        assert "SOG001\tG1A\tG1B\tG1C" in lines[1]

        # Check Row 2 (Paralog and Missing)
        # Sample_A should have G2A1,G2A2; Sample_C should be empty
        fields2 = lines[2].split("\t")
        assert fields2[0] == "SOG002"
        assert "G2A1" in fields2[1] and "G2A2" in fields2[1]
        assert "," in fields2[1]  # Comma separation for paralogs
        assert fields2[2] == "G2B"
        assert fields2[3] == ""  # Empty column for missing sample

    def test_write_empty_iterator(self, tmp_path):
        """Ensure the function writes at least a header if no SOGs are found."""
        output_file = tmp_path / "empty.tsv"
        write_og_table(iter([]), ["Sample_A"], output_file)

        lines = output_file.read_text().splitlines()
        assert len(lines) == 1
        assert lines[0] == "SOG_ID\tSample_A"
