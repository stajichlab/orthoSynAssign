# orthosynassign-vis

`orthosynassign-vis` is a companion visualization script used to verify the SOGs produced by [orthosynassign](orthosynassign.md). It uses [pyGenomeViz](https://github.com/moshi4/pyGenomeViz) to plot orthogroups and their synteny relationships, taking the original `orthogroup.tsv` alongside the SOG output and plotting selected SOGs using their previous orthogroup IDs as gene labels.

---

## Command Syntax

```sh
orthosynassign-vis --og_file OG_FILE --sog_file SOG_FILE --bed FILE [FILE ...] --sog SOG [SOG ...] [OPTIONS]
```

Example:

```sh
orthosynassign-vis --og_file orthogroups.tsv --sog_file Refined_SOGs.tsv --bed *.bed --sog SOG000023.OG0000025 -f svg
```

---

## Arguments & Options Reference

### Required Arguments

| Flag | Long Flag | Description |
| :--- | :--- | :--- |
| | `--og_file` | Path to the original `orthogroups.tsv` file. |
| | `--sog_file` | Path to the SOG (`Refined_SOGs.tsv`) file. |
| | `--bed` | Path of BED formatted genome annotation files. Space-separated list. |
| | `--sog` | Plot the SOG(s) produced by the previous `orthosynassign` analysis. Space-separated list of one or more SOG IDs. |

### Options

| Flag | Long Flag | Default | Description |
| :--- | :--- | :--- | :--- |
| `-w` | `--window` | `8` | The window size applied to the previous `orthosynassign` analysis. |
| `-o` | `--output` | `visualize_[sog_file]` | Output directory. |
| `-f` | `--fmt` | `png` | Output image format (`png`, `jpg`, `svg`, or `pdf`). |
| `-k` | `--keep_all_genes` | `False` | Keep genes that are not assigned to any orthogroup. |
| | `--collapse` | `tandem` | Collapse focal genes of the same genome into a single track (`tandem`, `window` or `none`). |
| `-c` | `--combine` | `False` | Combine all SOGs in a single figure file `combined.[fmt]`. |
| `-v` | `--verbose` | `False` | Enable verbose logging. |
| `-V` | `--version` | | Show program's version number and exit. |
| `-h` | `--help` | | Show help message and exit. |

:::{note}
`--window` should match the `-w/--window` value used in the original `orthosynassign` run being visualized, so that the plotted genomic context matches the window that was actually used for synteny support during partitioning.
:::

### Contig ends

A window stops either because `--window` is used up or because the contig ends, and both look the same in a figure. A dark bar is therefore drawn beside the first or the last gene of a track when that gene is the first or the last gene of its contig, so a short flank can be told from a truncated one. Circular genomes have no end and never get the mark.

---

### Collapsing focal genes

A SOG can hold several focal genes from the same genome, for instance a run of tandem copies. Each of them is drawn as its own track by default, which repeats nearly the same genomic context many times. `--collapse` merges such tracks:

| Mode | Behaviour |
| :--- | :--- |
| `tandem` (default) | Focal genes of the same genome that are drawn next to each other on the same contig are drawn as one track. Adjacency is judged on the genes actually plotted, so copies separated only by genes that the figure leaves out (genes without an orthogroup when `-k/--keep_all_genes` is off) still count as tandem. |
| `window` | Focal genes of the same genome whose plotted windows share at least one gene are drawn as one track. This is broader than `tandem` and also merges copies separated by a few genes. |

A track is only merged into another when the kept window actually shows its focal gene, so no focal gene disappears from the figure. A long chain of overlapping windows is therefore split into as few tracks as are needed to show every copy.
| `none` | No collapsing: one track per focal gene, as in v1.4.0 and earlier. |

Collapsing is transitive, so a run of three tandem copies becomes a single track. A collapsed track shows the union of the windows of the focal genes it stands for, so no gene is lost by merging, its title lists the collapsed gene IDs (truncated to three plus a count), and every focal gene of the group is drawn in the focal gene color `#fcfc42`.

```sh
orthosynassign-vis --og_file orthogroups.tsv --sog_file Refined_SOGs.tsv --bed *.bed \
  --sog SOG000023.OG0000025 --collapse window
```

---

`--sog` accepts one or more SOG IDs, each rendered as its own figure by default (`<output>/<sog_id>.<fmt>`). Pass `-c/--combine` to instead render all requested SOGs together into a single `<output>/combined.<fmt>` figure:

```sh
orthosynassign-vis --og_file orthogroups.tsv --sog_file Refined_SOGs.tsv --bed *.bed \
  --sog SOG000023.OG0000025 SOG000024.OG0000025 SOG000025.OG0000025 -c -f svg
```

---

## Examples

Plot multiple SOGs, each as its own figure:

```sh
orthosynassign-vis --og_file orthogroup.tsv --sog_file Refined_SOGs.tsv --bed *.bed --sog SOG000001.OG0000001 SOG000002.OG0000007
```

Combine multiple SOGs into a single figure:

```sh
orthosynassign-vis --og_file orthogroup.tsv --sog_file Refined_SOGs.tsv --bed *.bed --sog SOG000001.OG0000001 SOG000002.OG0000007 -c
```

Specify output directory for figures:

```sh
orthosynassign-vis --og_file orthogroup.tsv --sog_file Refined_SOGs.tsv --bed *.bed --sog SOG000001.OG0000001 -o fig
```

Specify window size applied in the previous orthosynassign analysis:

```sh
orthosynassign-vis --og_file orthogroup.tsv --sog_file Refined_SOGs.tsv --bed *.bed --sog SOG000001.OG0000001 -w 10
```

Save in svg format:

```sh
orthosynassign-vis --og_file orthogroup.tsv --sog_file Refined_SOGs.tsv --bed *.bed --sog SOG000001.OG0000001 -f svg
```

With verbose output:

```sh
orthosynassign-vis --og_file orthogroup.tsv --sog_file Refined_SOGs.tsv --bed *.bed --sog SOG000001.OG0000001 -v
```

---

## Reading the Output Figure

- Genes belonging to the observed SOG are labelled in **yellow**.
- Genes assigned to the same orthogroup within the given window are labelled in other chromatic colors.
- Genes with orthologs that are **not syntenic** within the given window are labelled in **gray**.
- Pass `-k/--keep_all_genes` to additionally display genes that were not assigned to any orthogroup.

### Single SOG Figures (default)

The image is written to `<output>/<sog_id>.<fmt>`, e.g. `visualize_Refined_SOGs/SOG000023.OG0000025.svg`:

```{image} ../../misc/SOG000023.OG0000025.svg
:alt: Syntenic orthogroup SOG000023
:width: 800px
```
```{image} ../../misc/SOG000024.OG0000025.svg
:alt: Syntenic orthogroup SOG000024
:width: 800px
```
```{image} ../../misc/SOG000025.OG0000025.svg
:alt: Syntenic orthogroup SOG000025
:width: 800px
```

### Combined Figure (`-c/--combine`)

The image is written to `<output>/combined.<fmt>`:

```{image} ../../misc/combined.svg
:alt: Combined figure of SOG000023, SOG000024, and SOG000025
:width: 800px
```
