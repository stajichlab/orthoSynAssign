# Preparing Input Files

orthoSynAssign requires two kinds of input: an **OrthoFinder orthogroup file** and one **BED-formatted genome annotation file per species**.

---

## 1. Orthogroup File

The `orthogroup.tsv` (or the `N0.tsv` file found under OrthoFinder's phylogenetic hierarchical orthogroups directory) must be tab-separated with:

- **First column**: Orthogroup ID (e.g., `OG0000001`)
- **Subsequent columns**: Protein IDs for each species, one column per species, with the column header being the species name

This is the standard `Orthogroups.tsv` produced by [OrthoFinder](https://github.com/davidemms/OrthoFinder).

---

## 2. BED Genome Annotation Files

Most genome annotations are distributed as GFF3 or GTF files. Because the 9th ("attribute") column of GFF3/GTF is highly flexible, it is often hard to reliably parse out the exact protein IDs used by OrthoFinder. orthoSynAssign instead consumes a standardized BED format that links genomic coordinates directly to protein IDs.

### Converting GFF3/GTF to BED with `gff2bed.py`

The utility script `misc/gff2bed.py` converts GFF3 (or GTF) files into this BED format:

```sh
python misc/gff2bed.py --collapse-transcripts genome.gff3 > genome.bed
```

```
usage: gff2bed.py [-h]
                  [--style {eupathdb,sequence_ontology,refseq,ensembl,gtf}]
                  [--format {auto,gff3,gtf}] [--list-styles] [--use-gene-id]
                  [--collapse-transcripts]
                  [gff_file]

Convert GFF3 or GTF to BED format, mapping proteins to their genes.

positional arguments:
  gff_file              GFF3 or GTF file (gzipped or plain text)

options:
  -h, --help            show this help message and exit
  --style {eupathdb,sequence_ontology,refseq,ensembl,gtf}
                        Format style (default: auto-selected from the detected dialect)
  --format {auto,gff3,gtf}
                        Input dialect (default: auto-detect from file content)
  --list-styles         List available format styles and exit
  --use-gene-id         Use gene ID as name column instead of transcript ID (funannotate only)
  --collapse-transcripts
                        Collapse all transcripts per gene into one line with semicolon-separated IDs (funannotate only)
```

Enable `--collapse-transcripts` to match orthoSynAssign's isoform-processing design: this collapses genes with multiple isoforms into a single BED entry, with the 4th column containing all associated protein IDs concatenated with a semicolon (`;`).

### Manual BED Preparation

:::{important}
If you prepare BED files manually (without `gff2bed.py`), you **must** separate multi-isoform protein IDs with a semicolon (`;`) in the 4th column. orthoSynAssign is specifically programmed to use this delimiter to resolve isoform-related mapping issues automatically.
:::

---

## Next Steps

Once you have your orthogroup file and one BED file per species, continue to [Quick Start](quickstart.md) or the [orthosynassign](usage/orthosynassign.md) command reference to run the partitioning.
