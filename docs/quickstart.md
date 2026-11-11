# Quick Start

This walkthrough uses the example dataset bundled in the [`example`](https://github.com/stajichlab/orthoSynAssign/tree/main/example) directory of the repository, which contains three fungal BED annotations and an orthogroup file:

```
FungiDB-68_AfumigatusA1163.bed
FungiDB-68_AfumigatusAf293.bed
FungiDB-68_AnovofumigatusIBT16806.bed
orthogroups.tsv
```

---

## 1. Run the Refinement

From within the `example` directory, run:

```sh
orthosynassign --og_file orthogroups.tsv --bed *.bed -o Refined_SOGs.tsv
```

This refines the OrthoFinder-style `orthogroups.tsv` using synteny information from the three BED files and writes the result to `Refined_SOGs.tsv`. See [orthosynassign](usage/orthosynassign.md) for all available options (window size, ratio threshold, thread count, etc.).

---

## 2. Visualize a Refined Orthogroup

To verify a specific refined syntenic orthogroup (SOG), e.g. `SOG000023.OG0000025`, plot it against its original OrthoFinder orthogroup context:

```sh
orthosynassign-vis --og_file orthogroups.tsv --sog_file Refined_SOGs.tsv --bed *.bed --sog SOG000023.OG0000025 -f svg
```

The figure is written to `visualize_Refined_SOGs/SOG000023.OG0000025.svg`. See [orthosynassign-vis](usage/orthosynassign-vis.md) for all available options, including `-c/--combine` to plot multiple SOGs in a single figure.

In the resulting figure:

- Genes belonging to the observed refined orthogroup are labelled in **yellow**.
- Genes assigned to the same orthogroup within the given window are labelled in other chromatic colors.
- Genes with orthologs that are **not syntenic** within the given window are labelled in **gray**.

```{image} ../misc/SOG000023.OG0000025.svg
:alt: A refined orthogroup SOG000023
:width: 800px
```

---

## Next Steps

- Read [Preparing Input Files](input-files.md) to prepare your own BED files from GFF3/GTF annotations.
- See [Pipeline Architecture & Benchmarks](architecture.md) to understand how the synteny-refinement algorithm works and how it scales.
