# Pipeline Architecture & Benchmarks

This page summarizes the methodology and benchmark results reported in the orthoSynAssign preprint, [orthoSynAssign: refine orthogroups using synteny information](https://www.biorxiv.org/content/10.64898/2026.08.10.744007v1) (Tsai, Pina Paez & Stajich, bioRxiv 2026).

---

## Motivation

Clustering-based orthogroup inference (e.g. [OrthoFinder]) often produces complex, many-to-many orthogroups that mix true one-to-one orthologs with confounding paralogs. Incorporating **synteny** — the conserved local gene-order context around a gene — provides a robust way to split these clusters into high-granularity, single-copy syntenic orthogroups (SOGs). orthoSynAssign is a from-scratch rewrite of [OrthoRefine], which pioneered this approach in C++ but suffered from sample-matching bugs and high memory usage at scale.

---

## Hybrid Python/Rust Architecture

orthoSynAssign pairs an accessible Python command-line interface — for parsing OrthoFinder orthogroup files, BED annotations, and generating CLI/visualization output — with a synteny-comparison engine implemented in [Rust]. This split keeps installation and scripting simple (Python) while the performance-critical, per-gene neighborhood comparisons run in compiled, memory-efficient Rust code operating on optimized genomic-range data structures.

At a high level, for every gene shared between two genomes in a candidate orthogroup, orthoSynAssign inspects the `-w/--window` neighboring genes on either side and checks how many of them also have a syntenic (co-located) counterpart in the other genome. If the fraction of supporting neighbors meets the `-r/--ratio_threshold`, the gene pair is retained as syntenic, refining the original orthogroup into one or more smaller SOGs. See [orthosynassign](usage/orthosynassign.md) for the corresponding CLI options.

---

## Multithreading Scalability Benchmark

To assess performance at scale, orthoSynAssign was benchmarked across **193 Eurotiomycetes genomes** (a fungal class including *Aspergillus* and related genera) while varying the thread count (`-t/--threads`):

- **Runtime** decreased substantially as thread count increased, confirming strong multithreaded scalability of the Rust synteny engine.
- **Memory usage** remained strictly bounded and effectively **thread-independent** — i.e., running with more threads did not multiply peak memory usage, an important property for deployment on shared HPC systems with per-job memory limits.

These results indicate orthoSynAssign is well suited to high-throughput phylogenomic workflows involving hundreds of genomes.

---

## See Also

- [Citation and References](citation.md)
- [orthosynassign](usage/orthosynassign.md) and [orthosynassign-vis](usage/orthosynassign-vis.md) command references

[OrthoFinder]: https://github.com/davidemms/OrthoFinder
[OrthoRefine]: https://github.com/jl02142/OrthoRefine
[Rust]: https://rust-lang.org/
