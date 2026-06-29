from __future__ import annotations

import logging
import sys
from argparse import Namespace
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from orthosynassign.lib.parsers import AnnotationParser
from orthosynassign.visualize import AUTHOR, VERSION, main, run_cli


@pytest.fixture
def args_factory(tmp_path) -> Namespace:
    """Returns a function that creates a mock args object with defaults."""

    def create_args(**kwargs):
        # Define your standard defaults here
        defaults = {
            "og_file": tmp_path / "default_og.tsv",
            "sog_file": tmp_path / "sog.tsv",
            "bed": tmp_path / "default.bed",
            "output": tmp_path / "output.tsv",
            "sog": ["SOG001", "SOG002"],
            "fmt": "svg",
            "keep_all_genes": False,
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


class TestMain:
    def test_returns_0_on_success(self, monkeypatch, args_factory, tmp_path):
        args = args_factory()

        monkeypatch.setattr("orthosynassign.visualize.setup_logging", lambda verbose: None)
        monkeypatch.setattr("orthosynassign.visualize.validate_annotations", lambda args: [])
        monkeypatch.setattr("orthosynassign.visualize.validate_orthogroup", lambda f: f)
        monkeypatch.setattr("orthosynassign.visualize.read_og_table", lambda f, g: [])
        monkeypatch.setattr("orthosynassign.visualize.get_visualize_engine", lambda g, o, s: MagicMock())

        result = main(args)
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

    def test_logs_start(self, monkeypatch, args_factory, caplog):
        args = args_factory()

        monkeypatch.setattr("orthosynassign.visualize.setup_logging", lambda verbose: None)
        monkeypatch.setattr("orthosynassign.visualize.validate_annotations", lambda args: [])
        monkeypatch.setattr("orthosynassign.visualize.validate_orthogroup", lambda f: f)
        monkeypatch.setattr("orthosynassign.visualize.read_og_table", lambda f, g: [])
        monkeypatch.setattr("orthosynassign.visualize.get_visualize_engine", lambda g, o, s: MagicMock())

        with caplog.at_level(logging.INFO):
            main(args)

        assert any("Starting Visualize" in record.message for record in caplog.records)

    def test_warns_on_missing_sog(self, monkeypatch, args_factory, caplog):
        args = args_factory(sog=["SOG_MISSING"])

        monkeypatch.setattr("orthosynassign.visualize.setup_logging", lambda verbose: None)
        monkeypatch.setattr("orthosynassign.visualize.validate_annotations", lambda args: [])
        monkeypatch.setattr("orthosynassign.visualize.validate_orthogroup", lambda f: f)
        monkeypatch.setattr("orthosynassign.visualize.read_og_table", lambda f, g: [])
        monkeypatch.setattr("orthosynassign.visualize.get_visualize_engine", lambda g, o, s: MagicMock())

        with caplog.at_level(logging.WARNING):
            result = main(args)

        assert result == 0
        assert any("SOG_MISSING" in record.message for record in caplog.records)

    def test_output_dir_created(self, monkeypatch, args_factory, tmp_path):
        output_dir = tmp_path / "my_output"
        args = args_factory(output=output_dir, sog=[])

        monkeypatch.setattr("orthosynassign.visualize.setup_logging", lambda verbose: None)
        monkeypatch.setattr("orthosynassign.visualize.validate_annotations", lambda args: [])
        monkeypatch.setattr("orthosynassign.visualize.validate_orthogroup", lambda f: f)
        monkeypatch.setattr("orthosynassign.visualize.read_og_table", lambda f, g: [])
        monkeypatch.setattr("orthosynassign.visualize.get_visualize_engine", lambda g, o, s: MagicMock())

        result = main(args)

        assert result == 0
        assert output_dir.exists()

    def test_default_output_dir_created(self, monkeypatch, args_factory, tmp_path):
        sog_file = tmp_path / "my_sogs.tsv"
        args = args_factory(output=None, sog_file=sog_file, sog=[])

        monkeypatch.setattr("orthosynassign.visualize.setup_logging", lambda verbose: None)
        monkeypatch.setattr("orthosynassign.visualize.validate_annotations", lambda args: [])
        monkeypatch.setattr("orthosynassign.visualize.validate_orthogroup", lambda f: f)
        monkeypatch.setattr("orthosynassign.visualize.read_og_table", lambda f, g: [])
        monkeypatch.setattr("orthosynassign.visualize.get_visualize_engine", lambda g, o, s: MagicMock())

        original_cwd = Path.cwd()
        import os

        os.chdir(tmp_path)
        try:
            result = main(args)
        finally:
            os.chdir(original_cwd)

        assert result == 0
        assert (tmp_path / f"visualize_{sog_file.stem}").exists()
