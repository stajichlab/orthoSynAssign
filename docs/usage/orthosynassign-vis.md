# orthosynassign-vis

`orthosynassign-vis` is a companion visualization script used to verify the refined results produced by [orthosynassign](orthosynassign.md). It uses [pyGenomeViz](https://github.com/moshi4/pyGenomeViz) to plot orthogroups and their synteny relationships, taking the original unrefined `orthogroup.tsv` alongside the refined orthogroup output and plotting selected refined orthogroups (SOGs) using their previous orthogroup IDs as gene labels.

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
| | `--sog_file` | Path to the refined orthogroups (`Refined_SOGs.tsv`) file. |
| | `--bed` | Path of BED formatted genome annotation files. Space-separated list. |
| | `--sog` | Plot the SOG(s) refined from the previous `orthosynassign` analysis. Space-separated list of one or more SOG IDs. |

### Options

| Flag | Long Flag | Default | Description |
| :--- | :--- | :--- | :--- |
| `-w` | `--window` | `8` | The window size applied to the previous `orthosynassign` analysis. |
| `-o` | `--output` | `visualize_[sog_file]` | Output directory. |
| `-f` | `--fmt` | `png` | Output image format (`png`, `jpg`, `svg`, or `pdf`). |
| `-k` | `--keep_all_genes` | `False` | Keep genes that are not assigned to any orthogroup. |
| `-c` | `--combine` | `False` | Combine all SOGs in a single figure file `combined.[fmt]`. |
| `-v` | `--verbose` | `False` | Enable verbose logging. |
| `-V` | `--version` | | Show program's version number and exit. |
| `-h` | `--help` | | Show help message and exit. |

:::{note}
`--window` should match the `-w/--window` value used in the original `orthosynassign` run being visualized, so that the plotted genomic context matches the window that was actually used for synteny support during refinement.
:::

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

- Genes belonging to the observed refined orthogroup (SOG) are labelled in **yellow**.
- Genes assigned to the same orthogroup within the given window are labelled in other chromatic colors.
- Genes with orthologs that are **not syntenic** within the given window are labelled in **gray**.
- Pass `-k/--keep_all_genes` to additionally display genes that were not assigned to any orthogroup.

### Single SOG Figures (default)

The image is written to `<output>/<sog_id>.<fmt>`, e.g. `visualize_Refined_SOGs/SOG000023.OG0000025.svg`:

```{image} ../../misc/SOG000023.OG0000025.svg
:alt: A refined orthogroup SOG000023
:width: 800px
```
```{image} ../../misc/SOG000024.OG0000025.svg
:alt: A refined orthogroup SOG000024
:width: 800px
```
```{image} ../../misc/SOG000025.OG0000025.svg
:alt: A refined orthogroup SOG000025
:width: 800px
```

### Combined Figure (`-c/--combine`)

The image is written to `<output>/combined.<fmt>`:

```{image} ../../misc/combined.svg
:alt: Combined figure of SOG000023, SOG000024, and SOG000025
:width: 800px
```
