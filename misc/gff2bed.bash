#!/bin/bash

# Check if a file was provided
if [ -z "$1" ]; then
    echo "Usage: $0 <annotation.gff.gz>"
    exit 1
fi

# For NCBI-formatted gff3
# entry_id="ID"
# gene_entry="gene"
# protein_entry="CDS"
# protein_id="protein_id"

# For example gff3
entry_id="ID"
gene_entry="protein_coding_gene"
protein_entry="CDS"
protein_id="protein_source_id"

awk -v entry_id=$entry_id -v gene_entry=$gene_entry -v protein_entry=$protein_entry -v protein_id=$protein_id -F'\t' '
    function get_attr(str, key) {
        split(str, attrs, ";")
        for (i in attrs) {
            if (attrs[i] ~ "^" key "=") {
                sub("^" key "=", "", attrs[i])
                return attrs[i]
            }
        }
        return ""
    }

    # PASS 1: Map everything
    NR==FNR {
        id = get_attr($9, entry_id)
        parent = get_attr($9, "Parent")

        if (id != "") {
            # Store the coordinates and scaffold for IDs (especially for genes)
            if ($3 == gene_entry) {
                gene_coords[id] = $1 "\t" ($4 - 1) "\t" $5
            }
            # Map this feature to its parent
            parent_map[id] = parent
        }

        # If it is a CDS, grab the protein_id and link it to its specific ID
        if ($3 == protein_entry) {
            pid = get_attr($9, protein_id)
            if (pid != "" && parent != "") {
                # We store the protein ID against the Parent (mRNA or Gene)
                # Using an associative array to keep unique protein IDs per parent
                raw_cds_to_pid[parent SUBSEP pid] = 1
            }
        }
        next
    }

    # PASS 2: Trace the lineage
    # We only care about entries that have a protein_id (stored in pass 1)
    END {
        for (pair in raw_cds_to_pid) {
            split(pair, parts, SUBSEP)
            current = parts[1]
            protein = parts[2]

            # Climb the tree: CDS_Parent -> mRNA -> Gene
            while (current != "" && !(current in gene_coords)) {
                current = parent_map[current]
            }

            # If we found the root Gene
            if (current in gene_coords) {
                # Check if this specific Protein was already added to this specific Gene
                # Using a flat compound key to avoid the "scalar as array" error
                if (!( (current SUBSEP protein) in seen_pairs )) {
                    if (gene_to_proteins[current] == "") {
                        gene_to_proteins[current] = protein
                    } else {
                        gene_to_proteins[current] = gene_to_proteins[current] ";" protein
                    }
                    seen_pairs[current SUBSEP protein] = 1
                }
            }
        }

        # Output the final BED records
        for (g_id in gene_to_proteins) {
            print gene_coords[g_id] "\t" gene_to_proteins[g_id]
        }
    }
' <(zcat $1) <(zcat $1) | sort -k 1,1 -k2,2n
