# orthoSynAssign

[![CI](https://github.com/stajichlab/orthoSynAssign/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/stajichlab/orthoSynAssign/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/stajichlab/orthoSynAssign/graph/badge.svg?token=mNxGMxekfv)](https://codecov.io/gh/stajichlab/orthoSynAssign)
[![Conda](https://anaconda.org/bioconda/orthosynassign/badges/version.svg)](https://anaconda.org/bioconda/orthosynassign)
[![License](https://img.shields.io/github/license/stajichlab/orthoSynAssign?label=license)](https://github.com/stajichlab/orthoSynAssign/blob/main/LICENSE)
[![DOI](https://zenodo.org/badge/1140709253.svg)](https://doi.org/10.5281/zenodo.18762979)

**orthoSynAssign** (Ortholog Synteny Assignment Tool) is a Python tool to partition orthogroups into syntenic orthogroups (SOGs) of positional orthologs using synteny information
inferred from genome annotation files (BED converted from GFF3/GTF). It is a Python re-implementation of
[OrthoRefine], addressing sample-matching and memory-usage issues in the original C++ tool, and is designed to be more
efficient, easier to use, and more flexible for custom analyses.

Described in the preprint: [orthoSynAssign: refine orthogroups using synteny information (bioRxiv, 2026)](https://www.biorxiv.org/content/10.64898/2026.08.10.744007v1).

---

## Key Features

- **Synteny-based Partitioning**: Splits complex, many-to-many OrthoFinder orthogroups into high-granularity, single-copy syntenic
  orthogroups (SOGs) by comparing gene-neighborhood context across genomes.
- **Hybrid Python/Rust Architecture**: An intuitive Python command-line interface wraps a core synteny-comparison engine written
  in [Rust], giving straightforward installation, seamless data parsing, and exceptional computational efficiency.
- **Isoform-aware Parsing**: Use semicolon-delimited multi-isoform protein IDs in BED annotation files to resolve isoform-related
  mapping issues automatically.
- **Multithreaded & Scalable**: Benchmarked across 193 Eurotiomycetes genomes with strong multithreading scalability, drastically
  reduced runtime, and a strictly bounded, thread-independent memory footprint.
- **Configurable Synteny Window**: Tunable window size and synteny-support ratio threshold let users control the stringency of
  synteny-based orthogroup splitting.
- **Companion Visualization Tool**: `orthosynassign-vis` uses [pyGenomeViz] to plot SOGs and their syntenic context
  for manual verification of results.

---

## Pipeline Flowchart

```{image} _static/images/flow_chart.svg
:alt: orthoSynAssign flowchart
:width: 800px
:align: center
```

---

## Pipeline Overview

1. Convert genome annotation files (GFF3/GTF) to BED format with `misc/gff2bed.py`, collapsing isoforms into semicolon-delimited
   protein ID lists.
2. Run `orthosynassign` on an OrthoFinder `Orthogroups.tsv` (or `N0.tsv`) file plus the BED annotation files to partition them into
   syntenic orthogroups (SOGs).
3. Optionally verify individual SOGs with `orthosynassign-vis`, which renders a synteny plot for a chosen SOG
   using [pyGenomeViz].

---

## Contents

- **[Installation](installation.md)** — Bioconda, Pixi, source, release tarball, and the developer setup.
- **[Preparing Input Files](input-files.md)** — OrthoFinder orthogroup files, converting GFF3/GTF to BED with `gff2bed.py`, and
  the isoform delimiter convention.
- **[Quick Start](quickstart.md)** — Walkthrough using the bundled `example` dataset.
- **Command reference**: [orthosynassign](usage/orthosynassign.md) · [orthosynassign-vis](usage/orthosynassign-vis.md)
- **[Pipeline Architecture & Benchmarks](architecture.md)** — Algorithm overview and multithreading scalability benchmarks from
  the preprint.
- **[Citation and References](citation.md)** — Publication metadata and citations for underlying dependencies.

```{toctree}
:caption: Getting started
:hidden:

installation
input-files
quickstart
```

```{toctree}
:caption: Command reference
:hidden:

usage/orthosynassign
usage/orthosynassign-vis
```

```{toctree}
:caption: Background
:hidden:

architecture
```

```{toctree}
:caption: Resources
:hidden:

citation
changelog
GitHub repository <https://github.com/stajichlab/orthoSynAssign>
Bioconda package <https://anaconda.org/bioconda/orthosynassign>
```

[OrthoRefine]: https://github.com/jl02142/OrthoRefine
[Rust]: https://rust-lang.org/
[pyGenomeViz]: https://github.com/moshi4/pyGenomeViz
