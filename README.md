# orthoSynAssign

Ortholog Synteny Assignment Tool - A Python tool for analyzing orthologous groups and synteny using OrthoFinder results and genome annotation files (GFF3/GTF).

## Features

- Parse GFF3 and GTF formatted genome annotation files
- Read OrthoFinder Orthogroups.tsv files
- Process multiple genome annotations from a directory
- Command-line interface for easy integration into pipelines
- Extensible library for custom analyses

## Installation

### From source

```bash
# Clone or navigate to the project directory
cd orthoSynAssign

# Install in development mode
pip install -e .

# Or install normally
pip install .
```

## Usage

### Command Line Interface

The main program `orthoSynAssign.py` provides a command-line interface:

```bash
# Using GFF3 files
python orthoSynAssign.py --gff_folder /path/to/gff_files --orthofinder Orthogroups.tsv

# Using GTF files
python orthoSynAssign.py --gtf_folder /path/to/gtf_files --orthofinder Orthogroups.tsv

# Specify output directory
python orthoSynAssign.py --gff_folder annotations/ --orthofinder Orthogroups.tsv -o results/

# Enable verbose logging
python orthoSynAssign.py --gff_folder annotations/ --orthofinder Orthogroups.tsv -v
```

### Command Line Options

- `--gff_folder PATH` - Path to folder containing GFF3 formatted genome annotation files
- `--gtf_folder PATH` - Path to folder containing GTF formatted genome annotation files  
  (Note: `--gff_folder` and `--gtf_folder` are mutually exclusive)
- `--orthofinder PATH` - Path to OrthoFinder Orthogroups.tsv file (required)
- `-o, --output PATH` - Output directory for results (default: output)
- `-v, --verbose` - Enable verbose logging

### Using as a Library

You can also use orthoSynAssign as a library in your own Python scripts:

```python
from orthoSynAssign.lib import read_gtf, read_gff3, read_orthofinder_table
from orthoSynAssign.lib.parsers import read_gff_folder, read_gtf_folder

# Read a single GFF3 file
features = read_gff3('genome.gff3')

# Read a single GTF file
features = read_gtf('genome.gtf')

# Read all GFF3 files from a folder
annotations = read_gff_folder('annotations/')

# Read OrthoFinder orthogroups
orthogroups = read_orthofinder_table('Orthogroups.tsv')

# Access the data
for species, feature_list in annotations.items():
    print(f"{species}: {len(feature_list)} features")

for og_id, species_genes in orthogroups.items():
    print(f"{og_id}: {len(species_genes)} species")
```

## File Formats

### GFF3/GTF Files

The tool expects standard GFF3 or GTF formatted genome annotation files with 9 tab-separated columns:
1. seqid - Chromosome/scaffold name
2. source - Source of the annotation
3. type - Feature type (gene, mRNA, CDS, etc.)
4. start - Start position (1-based)
5. end - End position (1-based, inclusive)
6. score - Score (or '.')
7. strand - Strand (+, -, or .)
8. phase - Phase (0, 1, 2, or .)
9. attributes - Semicolon-separated attributes

### OrthoFinder Orthogroups.tsv

The Orthogroups.tsv file from OrthoFinder should be tab-separated with:
- First column: Orthogroup ID (e.g., OG0000001)
- Subsequent columns: Gene IDs for each species (column headers are species names)

## Project Structure

```
orthoSynAssign/
├── orthoSynAssign/          # Main package
│   ├── __init__.py          # Package initialization
│   └── lib/                 # Library modules
│       ├── __init__.py      # Library initialization
│       └── parsers.py       # File parsing functions
├── orthoSynAssign.py        # Main command-line program
├── setup.py                 # Installation script
├── requirements.txt         # Python dependencies
├── .gitignore              # Git ignore file
└── README.md               # This file
```

## Development

### Adding New Features

The library is designed to be extensible. You can add new parsing functions or analysis methods to the `orthoSynAssign/lib/` directory.

### Testing

```bash
# Install development dependencies
pip install -e ".[dev]"

# Run tests (when implemented)
pytest tests/
```

## Requirements

- Python >= 3.8
- No external dependencies required for basic functionality

## License

[Add your license here]

## Citation

If you use orthoSynAssign in your research, please cite:

[Add citation information]

## Contact

[Add contact information]

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
