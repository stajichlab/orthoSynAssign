from __future__ import annotations

import logging
import sys
from argparse import Namespace
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from orthosynassign.visualize import AUTHOR, VERSION, _collapse_tracks, _track_label, main, run_cli


@pytest.fixture
def args_factory(tmp_path) -> Namespace:
    """Returns a function that creates a mock args object with defaults."""

    def create_args(**kwargs):
        # Define your standard defaults here
        defaults = {
            "og_file": tmp_path / "default_og.tsv",
            "sog_file": tmp_path / "sog.tsv",
            "bed": tmp_path / "default.bed",
            "sog": ["SOG001", "SOG002"],
            "window": 8,
            "output": tmp_path / "output.tsv",
            "fmt": "svg",
            "keep_all_genes": False,
            "collapse": "tandem",
            "combine": False,
            "verbose": True,
        }
        defaults.update(kwargs)
        return Namespace(**defaults)

    return create_args


# --- run_cli ---


class TestRunCli:
    @pytest.mark.parametrize("exit_code", [0, 1])
    def test_exit_code(self, monkeypatch: pytest.MonkeyPatch, exit_code: int):
        monkeypatch.setattr(
            sys, "argv", ["orthosynassign", "--og_file", "og", "--sog_file", "sog", "--bed", "bed", "--sog", "SOG001", "SOG002"]
        )
        monkeypatch.setattr("orthosynassign.visualize.main", lambda args: exit_code)

        with pytest.raises(SystemExit) as excinfo:
            run_cli()

        assert excinfo.value.code == exit_code

    def test_version(self, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture):
        # Force main to return exit code 0
        monkeypatch.setattr("orthosynassign.visualize.main", lambda args: 0)

        monkeypatch.setattr(sys, "argv", ["orthosynassign", "-V"])
        with pytest.raises(SystemExit) as excinfo:
            run_cli()

        out1 = capsys.readouterr().out
        assert excinfo.value.code == 0

        monkeypatch.setattr(sys, "argv", ["orthosynassign", "--version"])
        with pytest.raises(SystemExit) as excinfo:
            run_cli()

        out2 = capsys.readouterr().out

        assert excinfo.value.code == 0
        assert out1.strip() == out2.strip() == VERSION

    def test_help(self, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture):
        # Force main to return exit code 0
        monkeypatch.setattr("orthosynassign.visualize.main", lambda args: 0)

        monkeypatch.setattr(sys, "argv", ["orthosynassign", "-h"])
        with pytest.raises(SystemExit) as excinfo:
            run_cli()

        out1 = capsys.readouterr().out
        assert excinfo.value.code == 0

        monkeypatch.setattr(sys, "argv", ["orthosynassign", "--help"])
        with pytest.raises(SystemExit) as excinfo:
            run_cli()

        out2 = capsys.readouterr().out

        assert excinfo.value.code == 0
        assert out1 == out2
        assert f"Written by {AUTHOR}" in out1

    def test_missing_args(self, monkeypatch):
        monkeypatch.setattr(sys, "argv", ["orthosynassign", "--og_file", "og", "--sog_file", "sog", "--bed", "bed"])
        # Force main to return exit code 0
        monkeypatch.setattr("orthosynassign.visualize.main", lambda args: 0)

        with pytest.raises(SystemExit) as excinfo:
            run_cli()

        assert excinfo.value.code == 2

    def test_invalid_args(self, monkeypatch):
        monkeypatch.setattr(sys, "argv", ["orthosynassign", "--invalid", "og", "--bed", "bed"])
        # Force main to return exit code 0
        monkeypatch.setattr("orthosynassign.visualize.main", lambda args: 0)

        with pytest.raises(SystemExit) as excinfo:
            run_cli()

        assert excinfo.value.code == 2


# --- main ---
@pytest.fixture
def fake_main_env(monkeypatch, mock_genome_factory, mock_og_factory):
    monkeypatch.setattr("orthosynassign.visualize.setup_logging", lambda verbose: None)

    mock_anno1 = MagicMock()
    mock_anno1.parse.return_value = mock_genome_factory("G1")
    mock_anno2 = MagicMock()
    mock_anno2.parse.return_value = mock_genome_factory("G2")
    monkeypatch.setattr("orthosynassign.visualize.validate_annotations", lambda args: [mock_anno1, mock_anno2])

    monkeypatch.setattr("orthosynassign.visualize.validate_orthogroup", lambda f: f)
    monkeypatch.setattr(
        "orthosynassign.visualize.read_og_table", lambda f, g: [mock_og_factory("SOG001"), mock_og_factory("SOG002")]
    )

    mock_engine = MagicMock()
    mock_engine.get_aligned_og.return_value = []
    monkeypatch.setattr("orthosynassign.visualize.get_visualize_engine", lambda g, o, s: mock_engine)
    mock_render = MagicMock()
    monkeypatch.setattr("orthosynassign.visualize._render_sog_figure", mock_render)

    return {"render": mock_render}


class TestMain:
    def test_returns_0_on_success(self, fake_main_env, args_factory, tmp_path):
        args = args_factory()
        result = main(args)
        fake_main_env["render"].call_count == 2
        assert result == 0

    def test_returns_2_on_file_not_found(self, monkeypatch, args_factory):
        args = args_factory()

        monkeypatch.setattr("orthosynassign.visualize.setup_logging", lambda verbose: None)
        monkeypatch.setattr(
            "orthosynassign.visualize.validate_annotations",
            lambda args: (_ for _ in ()).throw(FileNotFoundError("File not found")),
        )

        result = main(args)
        assert result == 2

    def test_returns_1_on_generic_exception(self, monkeypatch, args_factory):
        args = args_factory()

        monkeypatch.setattr("orthosynassign.visualize.setup_logging", lambda verbose: None)
        monkeypatch.setattr(
            "orthosynassign.visualize.validate_annotations",
            lambda args: (_ for _ in ()).throw(Exception("Something went wrong")),
        )

        result = main(args)
        assert result == 1

    def test_returns_130_on_keyboard_interrupt(self, monkeypatch, args_factory):
        args = args_factory()

        monkeypatch.setattr("orthosynassign.visualize.setup_logging", lambda verbose: None)
        monkeypatch.setattr(
            "orthosynassign.visualize.validate_annotations",
            lambda args: (_ for _ in ()).throw(KeyboardInterrupt()),
        )

        result = main(args)
        assert result == 130

    def test_logs_start(self, fake_main_env, args_factory, caplog):
        args = args_factory()

        with caplog.at_level(logging.INFO):
            main(args)

        assert any("Starting Visualize" in record.message for record in caplog.records)

    def test_warns_on_missing_sog(self, fake_main_env, args_factory, caplog):
        args = args_factory(sog=["SOG_MISSING"])

        with caplog.at_level(logging.WARNING):
            result = main(args)

        assert result == 0
        assert any("SOG_MISSING" in record.message for record in caplog.records)

    def test_output_dir_created(self, fake_main_env, args_factory, tmp_path):
        output_dir = tmp_path / "my_output"
        args = args_factory(output=output_dir, sog=[])

        result = main(args)

        assert result == 0
        assert output_dir.exists()

    def test_default_output_dir_created(self, fake_main_env, args_factory, tmp_path):
        sog_file = tmp_path / "my_sogs.tsv"
        args = args_factory(output=None, sog_file=sog_file, sog=[])

        original_cwd = Path.cwd()
        import os

        os.chdir(tmp_path)
        try:
            result = main(args)
        finally:
            os.chdir(original_cwd)

        assert result == 0
        assert (tmp_path / f"visualize_{sog_file.stem}").exists()

    def test_combine(self, fake_main_env, args_factory, tmp_path):
        args = args_factory(combine=True)

        result = main(args)

        assert result == 0
        fake_main_env["render"].assert_called_once()


# --- _collapse_tracks ---


@pytest.fixture
def collapse_env(gene_factory, genome_factory):
    """Builds a genome whose genes 0..5 sit on one contig, plus a second genome.

    Returns:
        tuple: (genome_a, genome_b) with 6 and 2 genes respectively.
    """
    genome_a = genome_factory("A")
    for i in range(6):
        genome_a.add_gene(gene_factory(f"A{i}", "scaf1"))
    genome_b = genome_factory("B")
    for i in range(2):
        genome_b.add_gene(gene_factory(f"B{i}", "scaf1"))
    return genome_a, genome_b


class TestCollapseTracks:
    def test_none_keeps_every_track(self, collapse_env):
        genome_a, _ = collapse_env
        aligned = [(genome_a[2], [genome_a[1], genome_a[2]]), (genome_a[3], [genome_a[2], genome_a[3]])]

        result = _collapse_tracks(aligned, "none")

        assert [focal for focal, _, _ in result] == [genome_a[2], genome_a[3]]
        assert [members for _, _, members in result] == [[genome_a[2]], [genome_a[3]]]

    def test_tandem_collapses_adjacent_focal_genes(self, collapse_env):
        genome_a, _ = collapse_env
        # A2, A3 and A4 form one tandem run; A2 keeps the widest window
        aligned = [
            (genome_a[2], [genome_a[1], genome_a[2], genome_a[3], genome_a[4]]),
            (genome_a[3], [genome_a[2], genome_a[3]]),
            (genome_a[4], [genome_a[3], genome_a[4]]),
        ]

        result = _collapse_tracks(aligned, "tandem")

        assert len(result) == 1
        focal, window, members = result[0]
        assert focal is genome_a[2]
        assert len(window) == 4
        assert members == [genome_a[2], genome_a[3], genome_a[4]]

    def test_tandem_keeps_focal_genes_separated_by_a_plotted_gene_apart(self, collapse_env):
        genome_a, _ = collapse_env
        # A2 is plotted between A1 and A3, so they are not drawn next to each other
        aligned = [
            (genome_a[1], [genome_a[1], genome_a[2], genome_a[3]]),
            (genome_a[3], [genome_a[1], genome_a[2], genome_a[3]]),
        ]

        result = _collapse_tracks(aligned, "tandem")

        assert [focal for focal, _, _ in result] == [genome_a[1], genome_a[3]]

    def test_tandem_ignores_genes_left_out_of_the_figure(self, collapse_env):
        genome_a, _ = collapse_env
        # A2 exists in the genome but is not plotted (e.g. no orthogroup, -k off), so A1 and
        # A3 are drawn next to each other and must collapse
        aligned = [
            (genome_a[1], [genome_a[0], genome_a[1], genome_a[3]]),
            (genome_a[3], [genome_a[1], genome_a[3], genome_a[4]]),
        ]

        result = _collapse_tracks(aligned, "tandem")

        assert len(result) == 1
        assert result[0][2] == [genome_a[1], genome_a[3]]

    def test_tandem_keeps_focal_genes_on_different_contigs_apart(self, gene_factory, genome_factory):
        genome = genome_factory("A")
        genome.add_gene(gene_factory("A0", "scaf1"))
        genome.add_gene(gene_factory("A1", "scaf2"))
        aligned = [(genome[0], [genome[0]]), (genome[1], [genome[1]])]

        result = _collapse_tracks(aligned, "tandem")

        assert [focal for focal, _, _ in result] == [genome[0], genome[1]]

    def test_window_collapses_overlapping_windows(self, collapse_env):
        genome_a, _ = collapse_env
        # A1 and A3 are not adjacent, but their windows share A2
        aligned = [
            (genome_a[1], [genome_a[1], genome_a[2]]),
            (genome_a[3], [genome_a[2], genome_a[3]]),
        ]

        result = _collapse_tracks(aligned, "window")

        assert len(result) == 1
        assert result[0][2] == [genome_a[1], genome_a[3]]

    def test_window_keeps_disjoint_windows_apart(self, collapse_env):
        genome_a, _ = collapse_env
        aligned = [
            (genome_a[0], [genome_a[0], genome_a[1]]),
            (genome_a[4], [genome_a[4], genome_a[5]]),
        ]

        result = _collapse_tracks(aligned, "window")

        assert [focal for focal, _, _ in result] == [genome_a[0], genome_a[4]]

    def test_never_collapses_across_genomes(self, collapse_env):
        genome_a, genome_b = collapse_env
        aligned = [(genome_a[0], [genome_a[0], genome_a[1]]), (genome_b[0], [genome_b[0], genome_b[1]])]

        for mode in ("tandem", "window"):
            result = _collapse_tracks(aligned, mode)
            assert [focal for focal, _, _ in result] == [genome_a[0], genome_b[0]]

    def test_unknown_mode_raises(self, collapse_env):
        genome_a, _ = collapse_env
        aligned = [(genome_a[0], [genome_a[0]]), (genome_a[1], [genome_a[1]])]

        with pytest.raises(ValueError, match="Unknown collapse mode"):
            _collapse_tracks(aligned, "nope")


# --- _track_label ---


class TestTrackLabel:
    def test_single_gene(self, collapse_env):
        genome_a, _ = collapse_env
        assert _track_label([genome_a[0]]) == "A0"

    def test_truncates_long_groups(self, collapse_env):
        genome_a, _ = collapse_env
        assert _track_label(list(genome_a[0:5])) == "A0\nA1\nA2\n(+2 more)"
