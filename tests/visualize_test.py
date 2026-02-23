from __future__ import annotations

import sys
from argparse import Namespace
from pathlib import Path

import pytest

from orthosynassign.visualize import AUTHOR, VERSION, main, run_cli


@pytest.fixture
def args_factory() -> Namespace:
    """Returns a function that creates a mock args object with defaults."""

    def create_args(**kwargs):
        # Define your standard defaults here
        defaults = {
            "og_file": "default_og.tsv",
            "sog_file": "sog.tsv",
            "bed": "default.bed",
            "output": "output.tsv",
            "sog": ["SOG001", "SOG002"],
            "fmt": "svg",
            "keep_all_genes": False,
            "verbose": False,
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
def mock_visualize_dependencies(monkeypatch: pytest.MonkeyPatch):
    """Fixture to mock all heavy dependencies in the visualize module."""
    # Mocking external file/validation calls
    monkeypatch.setattr("orthosynassign.visualize.setup_logging", lambda x: None)
    monkeypatch.setattr("orthosynassign.visualize.validate_annotations", lambda x: [])
    monkeypatch.setattr("orthosynassign.visualize.validate_orthogroup", lambda x: x)
    monkeypatch.setattr("orthosynassign.visualize.read_og_table", lambda x, y: {})
    monkeypatch.setattr("orthosynassign.visualize._prepare_sog_visualization_data", lambda a, b, c, cpus: iter([]))
    monkeypatch.setattr("orthosynassign.visualize._render_sog_figure", lambda a, b, c: None)

    # Mock Path.mkdir and Path.unlink to prevent actual disk changes
    monkeypatch.setattr(Path, "mkdir", lambda *args, **kwargs: None)
    monkeypatch.setattr(Path, "unlink", lambda self, missing_ok=True: None)
    monkeypatch.setattr(Path, "replace", lambda self, target: None)


class TestMain:
    def test_main_success(self, args_factory, mock_visualize_dependencies):
        """Test that main returns 0 on a successful run."""
        args = args_factory()
        result = main(args)
        assert result == 0

    def test_main_keyboard_interrupt(self, monkeypatch, args_factory, mock_visualize_dependencies):
        """Test that main catches Ctrl+C and returns 1."""
        args = args_factory()

        def mock_interrupt(args):
            raise KeyboardInterrupt

        monkeypatch.setattr("orthosynassign.visualize.validate_annotations", mock_interrupt)

        result = main(args)
        assert result == 130

    def test_main_validation_error(self, monkeypatch, args_factory, mock_visualize_dependencies):
        """Test that main catches general errors and returns 1."""
        args = args_factory()

        def mock_crash(*args):
            raise FileNotFoundError

        monkeypatch.setattr("orthosynassign.visualize.read_og_table", mock_crash)

        result = main(args)
        assert result == 2

    def test_main_general_exception(self, monkeypatch, args_factory, mock_visualize_dependencies):
        """Test that main catches general errors and returns 1."""
        args = args_factory()

        def mock_crash(*args):
            raise RuntimeError

        monkeypatch.setattr("orthosynassign.visualize.read_og_table", mock_crash)

        result = main(args)
        assert result == 1
