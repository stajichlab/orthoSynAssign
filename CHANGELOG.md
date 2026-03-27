# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic
Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Added the `representative` attribute to the `Gene` class and the `is_isoform` keyword argument to the `add_gene` method of the
  `Genome` class, allowing for the management of gene isoforms.

- Added the `synteny` module with `SyntenyEngine` class to accelerate the analysis with array computing.

### Changed

- `calculate_synteny_ratio` now takes numpy array as input instead of a list of str, improving performance by utilizing array computing.

- Use disjoint set union (DSU) instead of BFS search for finding clusters after refinement. (`lib.dsu` module)

### Removed

- `SOG`, attribute and `Refine` method from `lib.orthogroup.Orthogroup` class.

- `compare_gene_sets` functions from `lib.orthogroup` module.

### Fixed

- Resolved an issue identified in [#2](https://github.com/stajichlab/orthoSynAssign/issues/2), specifically a bug within the `compare_gene_sets` function in the `orthogroup` module.

- Fixed the visualization script that was previously unable to correctly label and color orthogroups.


## [1.0.0] - 2026-02-24

### Added

- Translate OrthoRefine logic in Python.
- Multiprocessing support for refinement steps.
- Companion script to visualize the refined orthologs.

[Unreleased]: https://github.com/stajichlab/orthoSynAssign/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/stajichlab/orthoSynAssign/releases/tag/v1.0.0
