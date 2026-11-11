# orthosynassign

`orthosynassign` is the main program for running the synteny-based orthogroup partitioning analysis. It takes an OrthoFinder-style `Orthogroups.tsv` file (or `N0.tsv` under the phylogenetic hierarchical orthogroups directory) and outputs syntenic orthogroups (SOGs) of positional orthologs, using the genome annotation BED files described in [Preparing Input Files](../input-files.md).

---

## Command Syntax

```sh
orthosynassign --og_file OG_FILE --bed FILE [FILE ...] [OPTIONS]
```

Example:

```sh
orthosynassign --og_file orthogroups.tsv --bed *.bed -o Refined_SOGs.tsv
```

---

## Arguments & Options Reference

### Required Arguments

| Flag | Long Flag | Description |
| :--- | :--- | :--- |
| | `--og_file` | Path to OrthoFinder `Orthogroups.tsv` file. |
| | `--bed` | Path of BED formatted genome annotation files. Space-separated list. |

### Options

| Flag | Long Flag | Default | Description |
| :--- | :--- | :--- | :--- |
| `-w` | `--window` | `8` | Controls how many total genes are considered when determining synteny for a single gene. |
| `-r` | `--ratio_threshold` | `0.5` | Controls how many genes within a window must provide synteny support to classify the genes being compared as syntenous. |
| `-o` | `--output` | `Refined_SOGs-[YYYYMMDD-HHMMSS].tsv` (UTC timestamp) | Output of results. |
| `-t` | `--threads` | `4` | Number of CPUs to use. |
| `-v` | `--verbose` | `False` | Enable verbose logging. |
| `-V` | `--version` | | Show program's version number and exit. |
| `-h` | `--help` | | Show help message and exit. |

---

## Examples

Specify bed files separately:

```sh
orthosynassign --og_file orthogroup.tsv --bed file1.bed file2.bed file3.bed
```

Specify all bed files in a directory and processed in parallel with 6 CPUs:

```sh
orthosynassign --og_file orthogroup.tsv --bed *.bed -t 6
```

Specify output file name for results:

```sh
orthosynassign --og_file orthogroup.tsv --bed *.bed -o Refined_SOGs.tsv
```

Specify window size and ratio threshold:

```sh
orthosynassign --og_file orthogroup.tsv --bed *.bed -w 10 -r 0.8
```

With verbose output:

```sh
orthosynassign --og_file orthogroup.tsv --bed *.bed -v
```

---

## Tuning Window and Ratio Threshold

- **`--window`** sets how many neighboring genes (total, across the comparison) are inspected when evaluating whether a gene pair is syntenic. A larger window tolerates more genomic rearrangement between the compared genomes at the cost of more permissive matching.
- **`--ratio_threshold`** sets the minimum fraction of genes within that window that must show synteny support before a gene pair is classified as syntenic. Raising this threshold makes the partitioning stricter (fewer, more confident SOGs); lowering it is more permissive.

The defaults (`window=8`, `ratio_threshold=0.5`) are recommended starting points — see [Pipeline Architecture & Benchmarks](../architecture.md) for details on how the synteny-based partitioning algorithm works.

---

## Performance

The core synteny-comparison engine is implemented in [Rust], exposed to the `--threads` option for multithreaded execution. Benchmarks across 193 Eurotiomycetes genomes show strong scalability with additional threads and a memory footprint that stays bounded regardless of thread count — see [Pipeline Architecture & Benchmarks](../architecture.md) for details.

[Rust]: https://rust-lang.org/
