use numpy::PyReadonlyArray1;
use pyo3::prelude::*;
use std::collections::{HashMap, HashSet};

#[pyclass]
#[pyo3(name = "SyntenyEngine")]
pub struct SyntenyEngine {
    og_inputs: Vec<Vec<i32>>,
    seq_inputs: Vec<Vec<i16>>,
    orthogroups: Vec<Vec<(usize, usize)>>,
    shared_og_matrix: Vec<Vec<Vec<i32>>>,
}

#[pymethods]
impl SyntenyEngine {
    #[new]
    #[pyo3(text_signature = "(num_genomes, num_orthogroups, og_inputs, seq_inputs)")]
    pub fn new(
        num_orthogroups: usize,
        og_inputs: Vec<Vec<i32>>,
        seq_inputs: Vec<Vec<i16>>,
    ) -> PyResult<Self> {
        let mut orthogroups = vec![Vec::new(); num_orthogroups];

        for (genome_idx, (og_vec, _)) in og_inputs.iter().zip(&seq_inputs).enumerate() {
            for (gene_idx, &og_int) in og_vec.iter().enumerate() {
                if og_int >= 0 {
                    let og_idx_usize = og_int as usize;
                    if og_idx_usize < num_orthogroups {
                        orthogroups[og_idx_usize].push((genome_idx, gene_idx));
                    }
                }
            }
        }

        let shared_og_matrix = Self::build_shared_matrix(&og_inputs);

        Ok(SyntenyEngine {
            og_inputs,
            seq_inputs,
            orthogroups,
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
    fn build_shared_matrix(og_inputs: &[Vec<i32>]) -> Vec<Vec<Vec<i32>>> {
        let num_genomes = og_inputs.len();
        let genome_sets: Vec<HashSet<i32>> = og_inputs
            .iter()
            .map(|arr| arr.iter().filter(|&&id| id != -1).cloned().collect())
            .collect();

        let mut matrix = vec![vec![Vec::new(); num_genomes]; num_genomes];
        for i in 0..num_genomes {
            for j in i..num_genomes {
                let mut intersection: Vec<i32> = if i == j {
                    genome_sets[i].iter().cloned().collect()
                } else {
                    genome_sets[i]
                        .intersection(&genome_sets[j])
                        .cloned()
                        .collect()
                };
                intersection.sort_unstable();
                matrix[i][j] = intersection.clone();
                matrix[j][i] = intersection;
            }
        }
        matrix
    }

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

        let clusters = self.cluster_genes(refined_pairs, genes);
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
            get_window_logic(
                &self.seq_inputs[s_idx],
                s_gene_idx,
                window_size,
                |i| shared_ogs.binary_search(&self.og_inputs[s_idx][i]).is_ok(),
                &mut idx_buffer,
            );
            let mut win_ogs: Vec<i32> = idx_buffer
                .iter()
                .map(|&i| self.og_inputs[s_idx][i])
                .collect();
            win_ogs.sort_unstable();
            secondary_data.push((s_gene_idx, win_ogs));
        }

        for &p_gene_idx in primary_genes {
            get_window_logic(
                &self.seq_inputs[p_idx],
                p_gene_idx,
                window_size,
                |i| shared_ogs.binary_search(&self.og_inputs[p_idx][i]).is_ok(),
                &mut idx_buffer,
            );
            if idx_buffer.is_empty() {
                continue;
            }

            p_win_buffer.clear();
            p_win_buffer.extend(idx_buffer.iter().map(|&i| self.og_inputs[p_idx][i]));
            p_win_buffer.sort_unstable();

            let mut best_candidate = None;
            let mut max_r = -1.0;

            for (s_gene_idx, s_ogs) in &secondary_data {
                // Call function to calculate the ratio
                let ratio = calculate_synteny_ratio_logic(&p_win_buffer, s_ogs);
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

    fn cluster_genes(
        &self,
        pairs: Vec<((usize, usize), (usize, usize))>,
        all_genes: &[(usize, usize)],
    ) -> Vec<Vec<(usize, usize)>> {
        let n = all_genes.len();
        let gene_to_id: HashMap<(usize, usize), usize> =
            all_genes.iter().enumerate().map(|(i, &c)| (c, i)).collect();

        // data[i] < 0 => Root, value is -(rank + 1)
        let mut dsu = vec![-1; n];

        for (u, v) in pairs {
            if let (Some(&u_id), Some(&v_id)) = (gene_to_id.get(&u), gene_to_id.get(&v)) {
                let root_u = self.find_dsu(&mut dsu, u_id);
                let root_v = self.find_dsu(&mut dsu, v_id);

                if root_u != root_v {
                    // dsu[root] is negative.
                    // If dsu[root_u] is -1 and dsu[root_v] is -2:
                    // -1 > -2 is true, but -2 is the deeper tree.
                    if dsu[root_u] > dsu[root_v] {
                        // root_v is deeper, attach u to v
                        dsu[root_u] = root_v as i32;
                    } else if dsu[root_u] < dsu[root_v] {
                        // root_u is deeper, attach v to u
                        dsu[root_v] = root_u as i32;
                    } else {
                        // Ranks are equal, attach u to v and increment v's rank
                        dsu[root_u] = root_v as i32;
                        dsu[root_v] -= 1; // Rank becomes more negative
                    }
                }
            }
        }

        // Grouping remains the same, but now the root IDs will be consistent
        let mut clusters: HashMap<usize, Vec<(usize, usize)>> = HashMap::with_capacity(n);
        for i in 0..n {
            let r = self.find_dsu(&mut dsu, i);
            clusters.entry(r).or_default().push(all_genes[i]);
        }

        let mut result: Vec<Vec<(usize, usize)>> = clusters.into_values().collect();

        // Crucial for matching Python output exactly
        for c in &mut result {
            c.sort_unstable();
        }
        result.sort_unstable_by(|a, b| a[0].cmp(&b[0]));

        result
    }

    fn find_dsu(&self, dsu: &mut [i32], mut i: usize) -> usize {
        let mut root = i;
        while dsu[root] >= 0 {
            root = dsu[root] as usize;
        }
        while dsu[i] >= 0 {
            let n = dsu[i] as usize;
            dsu[i] = root as i32;
            i = n;
        }
        root
    }
}

fn get_window_logic<F>(
    seq: &[i16],
    gene_idx: usize,
    window_size: usize,
    og_is_valid: F,
    buffer: &mut Vec<usize>,
) where
    F: Fn(usize) -> bool,
{
    buffer.clear();
    let half_win = window_size / 2;
    let focal_seqid = seq[gene_idx];

    // Look Left: Find up to half_win valid indices before gene_idx

    let mut left_indices = Vec::with_capacity(half_win);
    let (mut i, mut left_count) = (gene_idx, 0);
    while i > 0 && left_count < half_win {
        i -= 1;
        if seq[i] == focal_seqid && og_is_valid(i) {
            left_indices.push(i);
            left_count += 1;
        }
    }
    // Since we scanned backwards, reverse to keep ascending order
    for &idx in left_indices.iter().rev() {
        buffer.push(idx);
    }

    // Look Right: Find up to half_win valid indices after gene_idx
    let (mut j, mut right_count) = (gene_idx, 0);
    while j < seq.len() - 1 && right_count < half_win {
        j += 1;
        if seq[j] == focal_seqid && og_is_valid(j) {
            buffer.push(j);
            right_count += 1;
        }
    }
}

#[pyfunction]
#[pyo3(name = "get_window")]
#[pyo3(text_signature = "(og_masked_array, seq_array, gene_idx, window_size)")]
pub fn get_window_py(
    og_masked_array: PyReadonlyArray1<bool>,
    seq_array: PyReadonlyArray1<i16>,
    gene_idx: usize,
    window_size: usize,
) -> PyResult<Vec<usize>> {
    // Convert NumPy to Rust Slices
    let og_masked = og_masked_array.as_slice()?;
    let seq = seq_array.as_slice()?;

    // Call the pure Rust logic
    let mut result_vec = Vec::with_capacity(window_size);

    get_window_logic(
        seq,
        gene_idx,
        window_size,
        |i| og_masked[i],
        &mut result_vec,
    );

    Ok(result_vec)
}

fn calculate_synteny_ratio_logic(win_a: &[i32], win_b: &[i32]) -> f64 {
    let len_a = win_a.len();
    let len_b = win_b.len();

    if len_a == 0 || len_b == 0 {
        return 0.0;
    }

    let mut matches = 0;
    let (mut i, mut j) = (0, 0);

    // Two-pointer walk on pre-sorted slices
    while i < len_a && j < len_b {
        if win_a[i] == win_b[j] {
            matches += 1;
            i += 1;
            j += 1;
        } else if win_a[i] < win_b[j] {
            i += 1;
        } else {
            j += 1;
        }
    }

    matches as f64 / std::cmp::max(len_a, len_b) as f64
}

#[pyfunction]
#[pyo3(name = "calculate_synteny_ratio")]
#[pyo3(text_signature = "(win_a, win_b)")]
pub fn calculate_synteny_ratio_py(mut win_a: Vec<i32>, mut win_b: Vec<i32>) -> PyResult<f64> {
    // 1. Sort them here so we can use the slice-based logic
    win_a.sort_unstable();
    win_b.sort_unstable();

    // 2. Call the maintainable slice-based helper
    Ok(calculate_synteny_ratio_logic(&win_a, &win_b))
}

#[pymodule]
fn rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // This registers your function so it appears as rs.get_window()
    m.add_class::<SyntenyEngine>()?;
    m.add_function(wrap_pyfunction!(get_window_py, m)?)?;
    m.add_function(wrap_pyfunction!(calculate_synteny_ratio_py, m)?)?;
    Ok(())
}
