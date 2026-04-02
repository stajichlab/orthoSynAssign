# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic
Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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

[Unreleased]: https://github.com/stajichlab/orthoSynAssign/compare/v1.1.0...HEAD

[1.1.0]: https://github.com/stajichlab/orthoSynAssign/compare/v1.0.0...v1.1.0

[1.0.0]: https://github.com/stajichlab/orthoSynAssign/releases/tag/v1.0.0
