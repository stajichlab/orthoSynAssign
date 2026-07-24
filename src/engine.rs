use crate::utils::*;
use pyo3::prelude::*;
use std::collections::HashMap;

#[pyclass]
#[pyo3(name = "SyntenyEngine")]
pub struct SyntenyEngine {
    orthogroups: Vec<Vec<(usize, usize)>>,
    ogs_vec: Vec<Vec<i32>>,
    seqids_vec: Vec<Vec<i16>>,
    circular_genome_vec: Vec<bool>,
    shared_og_matrix: Vec<Vec<Vec<i32>>>,
}

#[pymethods]
impl SyntenyEngine {
    #[new]
    pub fn new(
        num_ogs: usize,
        ogs_all: Vec<Vec<i32>>,
        seqids_all: Vec<Vec<i16>>,
        is_circular_all: Vec<bool>,
    ) -> PyResult<Self> {
        let orthogroups = get_orthogroups_vec(num_ogs, &ogs_all);
        let shared_og_matrix = build_shared_matrix(&ogs_all);

        Ok(SyntenyEngine {
            orthogroups,
            ogs_vec: ogs_all,
            seqids_vec: seqids_all,
            circular_genome_vec: is_circular_all,
            shared_og_matrix,
        })
    }

    #[pyo3(text_signature = "(self, og_idx, window_size, ratio_threshold)")]
    pub fn refine(
        &self,
        py: Python<'_>,
        og_idx: usize,
        window_size: usize,
        ratio_threshold: f64,
    ) -> Vec<Vec<(usize, usize)>> {
        py.detach(|| self.refine_logic(og_idx, window_size, ratio_threshold))
    }
}

/// Internal Rust-only Logic
impl SyntenyEngine {
    fn refine_logic(
        &self,
        og_idx: usize,
        window_size: usize,
        ratio_threshold: f64,
    ) -> Vec<Vec<(usize, usize)>> {
        let genes = &self.orthogroups[og_idx];
        if genes.is_empty() {
            return Vec::new();
        }

        // Group and IMMEDIATELY sort
        let mut genes_by_genome: HashMap<usize, Vec<usize>> = HashMap::new();
        for &(genome_idx, gene_idx) in genes {
            genes_by_genome
                .entry(genome_idx)
                .or_default()
                .push(gene_idx);
        }

        let mut genome_indices: Vec<usize> = genes_by_genome.keys().cloned().collect();
        genome_indices.sort_unstable();

        // CRITICAL: Sort the internal gene lists
        for genome_idx in &genome_indices {
            if let Some(list) = genes_by_genome.get_mut(genome_idx) {
                list.sort_unstable();
            }
        }

        let mut refined_pairs = Vec::new();

        // Deterministic Pairwise Loop
        for i in 0..genome_indices.len() {
            for j in i + 1..genome_indices.len() {
                let idx_a = genome_indices[i];
                let idx_b = genome_indices[j];
                let genes_a = &genes_by_genome[&idx_a];
                let genes_b = &genes_by_genome[&idx_b];
                let (p_idx, s_idx, p_genes, s_genes) = if genes_a.len() <= genes_b.len() {
                    (idx_a, idx_b, genes_a, genes_b)
                } else {
                    (idx_b, idx_a, genes_b, genes_a)
                };

                let pairs = self.compare_gene_pairs(
                    p_idx,
                    s_idx,
                    p_genes,
                    s_genes,
                    window_size,
                    ratio_threshold,
                );
                refined_pairs.extend(pairs);
            }
        }

        // Final Sorting of refined_pairs to ensure DSU input is identical
        refined_pairs.sort_unstable();

        let clusters = cluster_genes(refined_pairs, genes);
        clusters
            .into_iter()
            .filter(|cluster| cluster.len() > 1)
            .collect()
    }

    fn compare_gene_pairs(
        &self,
        p_idx: usize,
        s_idx: usize,
        primary_genes: &[usize],
        secondary_genes: &[usize],
        window_size: usize,
        ratio_threshold: f64,
    ) -> Vec<((usize, usize), (usize, usize))> {
        let shared_ogs = &self.shared_og_matrix[p_idx][s_idx];
        let mut idx_buffer = Vec::with_capacity(window_size);
        let mut p_win_buffer = Vec::with_capacity(window_size);
        let mut refined_pairs = Vec::new();

        let mut secondary_data = Vec::with_capacity(secondary_genes.len());
        for &s_gene_idx in secondary_genes {
            get_window(
                &self.seqids_vec[s_idx],
                s_gene_idx,
                window_size,
                self.circular_genome_vec[s_idx],
                |idx| shared_ogs.binary_search(&self.ogs_vec[s_idx][idx]).is_ok(),
                |idx| self.ogs_vec[s_idx][idx],
                &mut idx_buffer,
            );
            let mut win_ogs: Vec<i32> =
                idx_buffer.iter().map(|&i| self.ogs_vec[s_idx][i]).collect();
            win_ogs.sort_unstable();
            secondary_data.push((s_gene_idx, win_ogs));
        }

        for &p_gene_idx in primary_genes {
            get_window(
                &self.seqids_vec[p_idx],
                p_gene_idx,
                window_size,
                self.circular_genome_vec[p_idx],
                |idx| shared_ogs.binary_search(&self.ogs_vec[p_idx][idx]).is_ok(),
                |idx| self.ogs_vec[p_idx][idx],
                &mut idx_buffer,
            );
            if idx_buffer.is_empty() {
                continue;
            }

            p_win_buffer.clear();
            p_win_buffer.extend(idx_buffer.iter().map(|&i| self.ogs_vec[p_idx][i]));
            p_win_buffer.sort_unstable();

            let mut best_candidate = None;
            let mut max_r = -1.0;

            for (s_gene_idx, s_ogs) in &secondary_data {
                // Call function to calculate the ratio
                let ratio = calculate_synteny_ratio(&p_win_buffer, s_ogs);
                if ratio >= ratio_threshold - 1e-9 && ratio > max_r + 1e-9 {
                    max_r = ratio;
                    best_candidate = Some(*s_gene_idx);
                }
            }
            if let Some(s_best) = best_candidate {
                refined_pairs.push(((p_idx, p_gene_idx), (s_idx, s_best)));
            }
        }
        refined_pairs
    }
}

#[pyclass]
pub struct VisualizeEngine {
    orthogroups: Vec<Vec<(usize, usize)>>,
    ogs_vec: Vec<Vec<i32>>,
    seqids_vec: Vec<Vec<i16>>,
    circular_genome_vec: Vec<bool>,
    shared_og_matrix: Vec<Vec<Vec<i32>>>,
}

#[pymethods]
impl VisualizeEngine {
    #[new]
    pub fn new(
        sogs: Vec<Vec<(usize, usize)>>,
        ogs_all: Vec<Vec<i32>>,
        seqids_all: Vec<Vec<i16>>,
        is_circular_all: Vec<bool>,
    ) -> PyResult<Self> {
        let shared_og_matrix = build_shared_matrix(&ogs_all);

        Ok(VisualizeEngine {
            orthogroups: sogs,
            ogs_vec: ogs_all,
            seqids_vec: seqids_all,
            circular_genome_vec: is_circular_all,
            shared_og_matrix,
        })
    }

    pub fn get_aligned_og(
        &self,
        sog_idx: usize,
        window_size: usize,
        keep_all_genes: Option<bool>,
    ) -> Vec<((usize, usize), Vec<Option<usize>>)> {
        let og_windows = self.get_og_windows(sog_idx, window_size, keep_all_genes);
        align_windows(og_windows)
    }
}

impl VisualizeEngine {
    fn get_og_windows(
        &self,
        sog_idx: usize,
        window_size: usize,
        keep_all_genes: Option<bool>,
    ) -> Vec<((usize, usize), Vec<usize>)> {
        let genes = &self.orthogroups[sog_idx];
        if genes.is_empty() {
            return Vec::new();
        }

        let keep_all = keep_all_genes.unwrap_or(false);
        // (window_start_offset, window_end_offset)
        let mut boundaries: Vec<(i32, i32)> = vec![(0, 0); genes.len()];
        let mut idx_buffer = Vec::with_capacity(window_size);

        for i in 0..genes.len() {
            for j in i + 1..genes.len() {
                let (genome_a_idx, gene_a_idx) = genes[i];
                let (genome_b_idx, gene_b_idx) = genes[j];
                let shared = &self.shared_og_matrix[genome_a_idx][genome_b_idx];

                self.expand_boundary(
                    genome_a_idx,
                    gene_a_idx,
                    window_size,
                    &mut boundaries[i],
                    shared,
                    &mut idx_buffer,
                );
                self.expand_boundary(
                    genome_b_idx,
                    gene_b_idx,
                    window_size,
                    &mut boundaries[j],
                    shared,
                    &mut idx_buffer,
                );
            }
        }

        let mut result = Vec::with_capacity(genes.len());
        for (idx, &(genome_idx, focal_gene_idx)) in genes.iter().enumerate() {
            let (first, last) = boundaries[idx];
            let mut window_genes = Vec::new();
            let n_genes = self.seqids_vec[genome_idx].len() as i32;
            let is_circular = self.circular_genome_vec[genome_idx];

            for offset in first..=last {
                let gene_idx_i32 = focal_gene_idx as i32 + offset;
                let gene_idx = if is_circular {
                    gene_idx_i32.rem_euclid(n_genes) as usize
                } else {
                    if gene_idx_i32 < 0 || gene_idx_i32 >= n_genes {
                        continue;
                    }
                    gene_idx_i32 as usize
                };

                let current_og = self.ogs_vec[genome_idx][gene_idx];
                if keep_all || current_og >= 0 || gene_idx == focal_gene_idx {
                    window_genes.push(gene_idx);
                }
            }
            result.push(((genome_idx, focal_gene_idx), window_genes))
        }
        result
    }

    fn expand_boundary(
        &self,
        genome_idx: usize,
        gene_idx: usize,
        window_size: usize,
        boundary: &mut (i32, i32),
        shared: &[i32],
        // Avoid reapted mem allocation for performance
        idx_buffer: &mut Vec<usize>,
    ) {
        let n_genes = self.seqids_vec[genome_idx].len();
        let is_circular = self.circular_genome_vec[genome_idx];

        // Assuming get_window utility is defined elsewhere
        get_window(
            &self.seqids_vec[genome_idx],
            gene_idx,
            window_size,
            is_circular,
            |idx| shared.binary_search(&self.ogs_vec[genome_idx][idx]).is_ok(),
            |idx| self.ogs_vec[genome_idx][idx],
            idx_buffer,
        );

        if !idx_buffer.is_empty() {
            for &found_idx in idx_buffer.iter() {
                let mut diff = found_idx as i32 - gene_idx as i32;
                if is_circular {
                    let half = n_genes as i32 / 2;
                    if diff > half {
                        diff -= n_genes as i32;
                    } else if diff < -half {
                        diff += n_genes as i32;
                    }
                }
                boundary.0 = boundary.0.min(diff);
                boundary.1 = boundary.1.max(diff);
            }
        }
    }
}
