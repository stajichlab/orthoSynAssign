# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic
Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Added the `representative` attribute to the `Gene` class and the `is_isoform` keyword argument to the `add_gene` method of the
  `Genome` class, allowing for the management of gene isoforms.

### Fixed

- Resolved an issue identified in [#2](https://github.com/stajichlab/orthoSynAssign/issues/2), specifically a bug within the `compare_gene_sets` function in the `orthogroup` module.

## [1.0.0] - 2026-02-24

### Added

- Translate OrthoRefine logic in Python.
- Multiprocessing support for refinement steps.
- Companion script to visualize the refined orthologs.

[Unreleased]: https://github.com/stajichlab/orthoSynAssign/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/stajichlab/orthoSynAssign/releases/tag/v1.0.0
