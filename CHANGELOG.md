# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic
Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- When several genes of the secondary genome tie for the best synteny ratio, every one of them now gets a pairwise edge.
  Previously only the first in gene order was kept, so the result depended on gene order and a run of equally
  syntenic copies could be dropped (e.g. in YGOB OG0000006, `Suva_1.97`, `Suva_1.98` and `Suva_1.100` tie with
  `Suva_1.92` at 0.5 against the *S. cerevisiae* copies but were excluded). SOG membership can therefore differ from
  v1.4.0 on the same input: on the YGOB dataset (20 genomes, `-w 8 -r 0.5`) 122 genes are added, none are removed and
  23 SOGs are merged. Runtime is unchanged, as every candidate was already scored.
- Regenerated `example/SOGs.tsv` with the new rule: 6 more genes are assigned, the SOG count is unchanged.
- With ties kept, `utils::rescue_tandem_paralogs` no longer changes the result for linear genomes, because copies in a
  tandem run share the same window and tie. It is kept as a safety net for circular genomes.

## [1.4.0] - 2026-11-11

### Added

- Added a Sphinx documentation site published on Read the Docs, built from `docs/` with the MyST parser and the Furo theme and
  configured by `.readthedocs.yaml`. The pages are the former GitHub wiki content: installation, preparing input files, quick
  start, the `orthosynassign` and `orthosynassign-vis` command references, pipeline architecture and citation. The CHANGELOG is
  included as a page, and the new `docs` extra installs the build dependencies.

### Changed

- Tandem in-paralogs are now kept in SOGs, so SOG membership can differ from v1.3.1 on the same input. Previously, because
  each gene keeps only its single best syntenic match, the extra copies of a tandem array received no pairwise edge and were
  discarded as singletons. The new `utils::rescue_tandem_paralogs` runs after clustering and attaches such a copy to the cluster
  of an adjacent run-mate that was already validated, where the run is the consecutive genes of the same orthogroup on the same
  contig. Only physically adjacent copies are rescued, so distant paralogs are never chained, and the rescue can only move genes
  that are already members of the orthogroup being partitioned. On the bundled `example` dataset this adds 16 genes across 16
  SOGs, with no gene or SOG removed.
- The README, documentation and `CITATION.cff` now describe orthoSynAssign as partitioning orthogroups into syntenic orthogroups
  (SOGs) of positional orthologs, replacing the "refine/refinement" wording.
- Renamed the bundled example output from `example/Refined_SOGs.tsv` to `example/SOGs.tsv`.

## [1.3.1] - 2026-08-07

### Added

- Added the `-c/--combine` option to `visualize.py` to generate a single figure with multiple SOGs.

### Fixed

- Fixed a formatting issue in `CITATION.cff`.

- Fixed a bug in the Rust implementation of `utils::get_window_linear`.

## [1.3.0] - 2026-08-03

### Changed

- Skip tandem repeated genes during `get_window` to avoid disturbance of synteny calculation.

- Clean up the BED and orthogroups parsing by silencing the verbose log output.

## [1.2.0] - 2026-06-30

### Added

- Added Rust-implemented `VisualizeEngine` and utilize unified logic previously employed in the `SyntenyEngine` for visualization tasks.

- Added the get\_window logic for circular genome.

### Changed

- Refactored the `chromosome_type` attribute in the `Genome` class to `is_circular` for improved code clarity and simplified logic.

### Removed

- Python version object-oriented functions/methods for refinement and visualization have been superseded by Rust-implemented engine.

## [1.1.0] - 2026-04-02

### Added

- Added the `representative` attribute to the `Gene` class and the `is_isoform` keyword argument to the `add_gene` method of the
  `Genome` class, allowing for the management of gene isoforms.

- Added the `synteny` module to initiate a `SyntenyEngine` object written in Rust to accelerate the analysis and minimize memory
  overhead.

### Changed

- `calculate_synteny_ratio` and `get_window` are now implemented in Rust.

- Use disjoint set union (DSU) instead of BFS search for finding clusters after refinement.

- Use `ThreadPool` instead of `Pool` to minimize memory overhead.

- Sort the final refined clusters to ensure consistent results across runs.

### Removed

- `SOG`, attribute and `Refine` method from `lib.orthogroup.Orthogroup` class.

- `compare_gene_sets` functions from `lib.orthogroup` module.

### Fixed

- Resolved an issue identified in [#2](https://github.com/stajichlab/orthoSynAssign/issues/2), specifically a bug within the `compare_gene_sets` function in the `orthogroup` module. (incorporated into lib.rs)

- Fixed the visualization script that was previously unable to correctly label and color orthogroups.

## [1.0.0] - 2026-02-24

### Added

- Translate OrthoRefine logic in Python.
- Multiprocessing support for refinement steps.
- Companion script to visualize the refined orthologs.

[Unreleased]: https://github.com/stajichlab/orthoSynAssign/compare/v1.4.0...HEAD

[1.4.0]: https://github.com/stajichlab/orthoSynAssign/compare/v1.3.1...v1.4.0

[1.3.1]: https://github.com/stajichlab/orthoSynAssign/compare/v1.3.0...v1.3.1

[1.3.0]: https://github.com/stajichlab/orthoSynAssign/compare/v1.2.0...v1.3.0

[1.2.0]: https://github.com/stajichlab/orthoSynAssign/compare/v1.1.0...v1.2.0

[1.1.0]: https://github.com/stajichlab/orthoSynAssign/compare/v1.0.0...v1.1.0

[1.0.0]: https://github.com/stajichlab/orthoSynAssign/releases/tag/v1.0.0
