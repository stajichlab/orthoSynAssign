"""
Tests for GFF3 and GTF parsers.

Tests reading and parsing GFF3 and GTF files in the tests/data/file_formats directory.
Verifies that the expected number of features are found and parsed correctly.
"""

import unittest
import os
from pathlib import Path
from orthoSynAssign.lib.parsers import read_gff3, read_gtf, GFFFeature


class TestGFFParsers(unittest.TestCase):
    """Test suite for GFF3 and GTF file parsing."""

    @classmethod
    def setUpClass(cls):
        """Set up test data paths."""
        cls.data_dir = Path(__file__).parent / "data" / "file_formats"

        # GFF files (note: .gff.gz is actually GTF format, not GFF3)

        cls.gtf_files = [
            "GCF_000865725.1_ViralMultiSegProj15521_genomic.gtf.gz",
            "Nempa1_GeneCatalog_genes_20130904.gff.gz",
            "Pacta1_2_GeneCatalog_genes_20110630.gff.gz",
            "Rozal_SC1_GeneCatalog_genes_20141014.gff.gz",
        ]

        # GFF3 files (explicit .gff3.gz)
        cls.gff3_files = [
            "GCF_000865725.1_ViralMultiSegProj15521_genomic.gff.gz",
            "Pacta1_2.filtered_proteins.Filteredmodels2.gff3.gz",
        ]

    def test_gff3_files_exist(self):
        """Verify that GFF3 test data files exist."""
        for gff3_file in self.gff3_files:
            file_path = self.data_dir / gff3_file
            self.assertTrue(file_path.exists(), f"GFF3 file not found: {file_path}")

    def test_gtf_files_exist(self):
        """Verify that GTF test data files exist."""
        for gtf_file in self.gtf_files:
            file_path = self.data_dir / gtf_file
            self.assertTrue(file_path.exists(), f"GTF file not found: {file_path}")

    def test_parse_gff3_file(self):
        """Test parsing of GFF3 files."""
        for gff3_file in self.gff3_files:
            file_path = self.data_dir / gff3_file

            with self.subTest(file=gff3_file):
                features = read_gff3(str(file_path))

                # Verify we got features
                self.assertGreater(
                    len(features), 0, f"No features parsed from {gff3_file}"
                )

                # Verify all features are GFFFeature objects
                for feature in features:
                    self.assertIsInstance(feature, GFFFeature)
                    self.assertIsNotNone(feature.seqid)
                    self.assertIsNotNone(feature.type)
                    self.assertGreaterEqual(feature.start, 1)
                    self.assertGreaterEqual(feature.end, feature.start)

                print(f"✓ {gff3_file}: {len(features)} features parsed")

    def test_parse_gtf_file(self):
        """Test parsing of GTF files."""
        for gtf_file in self.gtf_files:
            file_path = self.data_dir / gtf_file

            with self.subTest(file=gtf_file):
                features = read_gtf(str(file_path))

                # Verify we got features
                self.assertGreater(
                    len(features), 0, f"No features parsed from {gtf_file}"
                )

                # Verify all features are GFFFeature objects
                for feature in features:
                    self.assertIsInstance(feature, GFFFeature)
                    self.assertIsNotNone(feature.seqid)
                    self.assertIsNotNone(feature.type)
                    self.assertGreaterEqual(feature.start, 1)
                    self.assertGreaterEqual(feature.end, feature.start)

                print(f"✓ {gtf_file}: {len(features)} features parsed")

    def test_gff3_feature_count_summary(self):
        """Test and report feature counts for all GFF3 files."""
        print("\n=== GFF3 File Feature Counts ===")
        for gff3_file in self.gff3_files:
            file_path = self.data_dir / gff3_file
            if not file_path.exists():
                continue

            features = read_gff3(str(file_path))
            feature_types = {}
            for feature in features:
                feature_types[feature.type] = feature_types.get(feature.type, 0) + 1

            print(f"\n{gff3_file}:")
            print(f"  Total features: {len(features)}")
            for ftype, count in sorted(feature_types.items()):
                print(f"    {ftype}: {count}")

    def test_gtf_feature_count_summary(self):
        """Test and report feature counts for all GTF files."""
        print("\n=== GTF File Feature Counts ===")
        for gtf_file in self.gtf_files:
            file_path = self.data_dir / gtf_file
            if not file_path.exists():
                continue

            features = read_gtf(str(file_path))
            feature_types = {}
            for feature in features:
                feature_types[feature.type] = feature_types.get(feature.type, 0) + 1

            print(f"\n{gtf_file}:")
            print(f"  Total features: {len(features)}")
            for ftype, count in sorted(feature_types.items()):
                print(f"    {ftype}: {count}")

    def test_gtf_features_have_valid_attributes(self):
        """Test that parsed GTF features have properly parsed attributes."""
        for gtf_file in self.gtf_files:
            file_path = self.data_dir / gtf_file
            if not file_path.exists():
                continue

            features = read_gtf(str(file_path))

            with self.subTest(file=gtf_file):
                # GTF files typically have attributes like gene_id, transcript_id
                features_with_attrs = [f for f in features if f.attributes]
                self.assertGreater(
                    len(features_with_attrs),
                    0,
                    f"No features with attributes found in {gtf_file}",
                )

                # Verify attribute dictionaries are properly formed
                for feature in features_with_attrs[:5]:  # Check first 5
                    self.assertIsInstance(feature.attributes, dict)


if __name__ == "__main__":
    unittest.main(verbosity=2)
