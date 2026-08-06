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
                &self.ogs_vec[s_idx],
                s_gene_idx,
                shared_ogs,
                window_size,
                self.circular_genome_vec[s_idx],
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
                &self.ogs_vec[p_idx],
                p_gene_idx,
                shared_ogs,
                window_size,
                self.circular_genome_vec[p_idx],
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
        shared_ogs: &[i32],
        // Avoid reapted mem allocation for performance
        idx_buffer: &mut Vec<usize>,
    ) {
        let n_genes = self.seqids_vec[genome_idx].len();
        let is_circular = self.circular_genome_vec[genome_idx];

        // Assuming get_window utility is defined elsewhere
        get_window(
            &self.seqids_vec[genome_idx],
            &self.ogs_vec[genome_idx],
            gene_idx,
            shared_ogs,
            window_size,
            is_circular,
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

// ... existing code ...

#[cfg(test)]
mod tests {
    use super::*;
    use pyo3::Python;

    // --- Shared Test Fixtures ---

    /// Two genomes, 6 genes each, sharing OGs 0-3.
    /// Genome 0: [0, 1, 2, 3, -1, -1]
    /// Genome 1: [0, 1, 2, 3, -1, -1]
    fn two_genome_engine() -> SyntenyEngine {
        let ogs_all = vec![vec![0i32, 1, 2, 3, -1, -1], vec![0i32, 1, 2, 3, -1, -1]];
        let seqids_all = vec![vec![0i16; 6], vec![0i16; 6]];
        let is_circular_all = vec![false, false];
        SyntenyEngine::new(4, ogs_all, seqids_all, is_circular_all).unwrap()
    }

    /// Same data but shuffled in genome 1 to test ordering robustness.
    /// Genome 0: [0, 1, 2, 3, -1]
    /// Genome 1: [3, 2, 1, 0, -1]
    fn shuffled_genome_engine() -> SyntenyEngine {
        let ogs_all = vec![vec![0i32, 1, 2, 3, -1], vec![3i32, 2, 1, 0, -1]];
        let seqids_all = vec![vec![0i16; 5], vec![0i16; 5]];
        let is_circular_all = vec![false, false];
        SyntenyEngine::new(4, ogs_all, seqids_all, is_circular_all).unwrap()
    }

    /// OG 0 present in only one genome — no pairs possible.
    fn single_genome_engine() -> SyntenyEngine {
        let ogs_all = vec![vec![0i32, 1, 2], vec![-1i32, -1, -1]];
        let seqids_all = vec![vec![0i16; 3], vec![0i16; 3]];
        let is_circular_all = vec![false, false];
        SyntenyEngine::new(3, ogs_all, seqids_all, is_circular_all).unwrap()
    }

    fn simple_visualize_engine() -> VisualizeEngine {
        // SOG 0: gene 0 from each genome (both have OG 0)
        let sogs = vec![vec![(0usize, 0usize), (1usize, 0usize)]];
        let ogs_all = vec![vec![0i32, 1, 2, 3, -1], vec![0i32, 1, 2, 3, -1]];
        let seqids_all = vec![vec![0i16; 5], vec![0i16; 5]];
        let is_circular_all = vec![false, false];
        VisualizeEngine::new(sogs, ogs_all, seqids_all, is_circular_all).unwrap()
    }

    fn circular_visualize_engine() -> VisualizeEngine {
        let sogs = vec![vec![(0usize, 0usize), (1usize, 0usize)]];
        let ogs_all = vec![vec![0i32, 1, 2, 3], vec![0i32, 1, 2, 3]];
        let seqids_all = vec![vec![0i16; 4], vec![0i16; 4]];
        let is_circular_all = vec![true, true];
        VisualizeEngine::new(sogs, ogs_all, seqids_all, is_circular_all).unwrap()
    }

    // -----------------------------------------------------------------------
    // SyntenyEngine::new
    // -----------------------------------------------------------------------

    #[test]
    fn test_synteny_engine_new_succeeds() {
        let result = SyntenyEngine::new(
            4,
            vec![vec![0i32, 1, 2, 3], vec![0i32, 1, 2, 3]],
            vec![vec![0i16; 4], vec![0i16; 4]],
            vec![false, false],
        );
        assert!(result.is_ok());
    }

    #[test]
    fn test_synteny_engine_new_empty_inputs() {
        // Zero genomes, zero OGs — should not panic
        let result = SyntenyEngine::new(0, vec![], vec![], vec![]);
        assert!(result.is_ok());
    }

    #[test]
    fn test_synteny_engine_new_single_genome() {
        let result = SyntenyEngine::new(2, vec![vec![0i32, 1]], vec![vec![0i16; 2]], vec![false]);
        assert!(result.is_ok());
    }

    #[test]
    fn test_synteny_engine_new_no_shared_ogs() {
        // Genomes with completely disjoint OG sets
        let result = SyntenyEngine::new(
            4,
            vec![vec![0i32, 1], vec![2i32, 3]],
            vec![vec![0i16; 2], vec![0i16; 2]],
            vec![false, false],
        );
        assert!(result.is_ok());
    }

    // -----------------------------------------------------------------------
    // SyntenyEngine::refine
    // -----------------------------------------------------------------------

    #[test]
    fn test_refine_matching_genomes_produces_cluster() {
        Python::attach(|py| {
            let engine = two_genome_engine();
            // OG 0 appears at position 0 in both genomes with identical neighbours
            let clusters = engine.refine(py, 0, 3, 0.3);
            assert!(
                !clusters.is_empty(),
                "Expected at least one syntenic cluster for a well-conserved OG"
            );
        });
    }

    #[test]
    fn test_refine_cluster_contains_one_gene_per_genome() {
        Python::attach(|py| {
            let engine = two_genome_engine();
            let clusters = engine.refine(py, 0, 3, 0.3);
            for cluster in &clusters {
                // Every cluster member must come from a different genome
                let genomes: Vec<usize> = cluster.iter().map(|&(g, _)| g).collect();
                let unique: std::collections::HashSet<_> = genomes.iter().collect();
                assert_eq!(
                    unique.len(),
                    genomes.len(),
                    "Cluster should not contain two genes from the same genome: {:?}",
                    cluster
                );
            }
        });
    }

    #[test]
    fn test_refine_all_clusters_have_more_than_one_member() {
        Python::attach(|py| {
            let engine = two_genome_engine();
            let clusters = engine.refine(py, 0, 3, 0.3);
            for cluster in &clusters {
                assert!(
                    cluster.len() > 1,
                    "Refine must only return clusters with >1 gene, got: {:?}",
                    cluster
                );
            }
        });
    }

    #[test]
    fn test_refine_empty_og_returns_empty() {
        Python::attach(|py| {
            // num_ogs = 4 but we request og_idx 3 which has no genes in seqids
            let engine = SyntenyEngine::new(
                5,
                vec![vec![0i32, 1, 2, 3], vec![0i32, 1, 2, 3]],
                vec![vec![0i16; 4], vec![0i16; 4]],
                vec![false, false],
            )
            .unwrap();
            // OG index 4 has no genes assigned
            let clusters = engine.refine(py, 4, 3, 0.3);
            assert!(
                clusters.is_empty(),
                "An OG with no genes should yield no clusters"
            );
        });
    }

    #[test]
    fn test_refine_single_genome_og_returns_empty() {
        Python::attach(|py| {
            let engine = single_genome_engine();
            // OG 0 only exists in genome 0 — no pairs to form
            let clusters = engine.refine(py, 0, 3, 0.3);
            assert!(
                clusters.is_empty(),
                "OG present in only one genome should produce no clusters"
            );
        });
    }

    #[test]
    fn test_refine_high_threshold_reduces_clusters() {
        Python::attach(|py| {
            let engine = shuffled_genome_engine();
            let loose = engine.refine(py, 0, 3, 0.1);
            let strict = engine.refine(py, 0, 3, 1.0);
            assert!(
                strict.len() <= loose.len(),
                "Stricter threshold should yield fewer or equal clusters"
            );
        });
    }

    #[test]
    fn test_refine_is_deterministic() {
        Python::attach(|py| {
            let engine = two_genome_engine();
            let first = engine.refine(py, 1, 3, 0.3);
            let second = engine.refine(py, 1, 3, 0.3);
            assert_eq!(first, second, "Refine must be deterministic");
        });
    }

    #[test]
    fn test_refine_window_size_one_still_runs() {
        Python::attach(|py| {
            let engine = two_genome_engine();
            // window_size = 1 is an edge case — should not panic
            let clusters = engine.refine(py, 0, 1, 0.3);
            // Result may be empty or not; just assert no panic
            let _ = clusters;
        });
    }

    #[test]
    fn test_refine_circular_genome() {
        Python::attach(|py| {
            let engine = SyntenyEngine::new(
                4,
                vec![vec![0i32, 1, 2, 3], vec![0i32, 1, 2, 3]],
                vec![vec![0i16; 4], vec![0i16; 4]],
                vec![true, true], // circular
            )
            .unwrap();
            // Should handle wrap-around without panic
            let clusters = engine.refine(py, 0, 3, 0.3);
            let _ = clusters;
        });
    }

    // -----------------------------------------------------------------------
    // VisualizeEngine::new
    // -----------------------------------------------------------------------

    #[test]
    fn test_visualize_engine_new_succeeds() {
        let result = VisualizeEngine::new(
            vec![vec![(0, 0), (1, 0)]],
            vec![vec![0i32, 1, 2], vec![0i32, 1, 2]],
            vec![vec![0i16; 3], vec![0i16; 3]],
            vec![false, false],
        );
        assert!(result.is_ok());
    }

    #[test]
    fn test_visualize_engine_new_empty_sogs() {
        let result = VisualizeEngine::new(vec![], vec![], vec![], vec![]);
        assert!(result.is_ok());
    }

    #[test]
    fn test_visualize_engine_new_empty_sog_entry() {
        // A SOG list that contains an empty inner vec
        let result = VisualizeEngine::new(
            vec![vec![]],
            vec![vec![0i32, 1]],
            vec![vec![0i16; 2]],
            vec![false],
        );
        assert!(result.is_ok());
    }

    // -----------------------------------------------------------------------
    // VisualizeEngine::get_aligned_og
    // -----------------------------------------------------------------------

    #[test]
    fn test_get_aligned_og_empty_sog_returns_empty() {
        let engine = VisualizeEngine::new(
            vec![vec![]],
            vec![vec![0i32, 1]],
            vec![vec![0i16; 2]],
            vec![false],
        )
        .unwrap();
        let result = engine.get_aligned_og(0, 3, None);
        assert!(result.is_empty(), "Empty SOG should produce empty output");
    }

    #[test]
    fn test_get_aligned_og_result_length_matches_sog_size() {
        let engine = simple_visualize_engine();
        let sog_size = 2; // two genes in SOG 0
        let result = engine.get_aligned_og(0, 3, None);
        assert_eq!(
            result.len(),
            sog_size,
            "Output should have one entry per gene in the SOG"
        );
    }

    #[test]
    fn test_get_aligned_og_gene_identifiers_match_sog() {
        let engine = simple_visualize_engine();
        let result = engine.get_aligned_og(0, 3, None);
        let identifiers: Vec<(usize, usize)> = result.iter().map(|&(id, _)| id).collect();
        assert!(
            identifiers.contains(&(0, 0)),
            "Result should contain gene (genome=0, gene=0)"
        );
        assert!(
            identifiers.contains(&(1, 0)),
            "Result should contain gene (genome=1, gene=0)"
        );
    }

    #[test]
    fn test_get_aligned_og_windows_same_length() {
        let engine = simple_visualize_engine();
        let result = engine.get_aligned_og(0, 3, None);
        // align_windows should pad all rows to equal length
        let lengths: Vec<usize> = result.iter().map(|(_, w)| w.len()).collect();
        assert!(
            lengths.windows(2).all(|p| p[0] == p[1]),
            "All aligned windows must have the same length, got {:?}",
            lengths
        );
    }

    #[test]
    fn test_get_aligned_og_keep_all_genes_true_vs_false() {
        let engine = simple_visualize_engine();
        let with_all = engine.get_aligned_og(0, 3, Some(true));
        let og_only = engine.get_aligned_og(0, 3, Some(false));
        // keep_all=true may introduce extra None slots for non-OG genes;
        // windows should be at least as wide
        let max_with_all: usize = with_all.iter().map(|(_, w)| w.len()).max().unwrap_or(0);
        let max_og_only: usize = og_only.iter().map(|(_, w)| w.len()).max().unwrap_or(0);
        assert!(
            max_with_all >= max_og_only,
            "keep_all=true windows should be >= keep_all=false"
        );
    }

    #[test]
    fn test_get_aligned_og_none_keep_all_behaves_like_false() {
        let engine = simple_visualize_engine();
        let with_none = engine.get_aligned_og(0, 3, None);
        let with_false = engine.get_aligned_og(0, 3, Some(false));
        assert_eq!(
            with_none, with_false,
            "keep_all_genes=None should behave identically to Some(false)"
        );
    }

    #[test]
    fn test_get_aligned_og_circular_does_not_panic() {
        let engine = circular_visualize_engine();
        // Gene 0 with window_size=3 on a 4-gene circular genome should wrap
        let result = engine.get_aligned_og(0, 3, None);
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_get_aligned_og_is_deterministic() {
        let engine = simple_visualize_engine();
        let first = engine.get_aligned_og(0, 3, None);
        let second = engine.get_aligned_og(0, 3, None);
        assert_eq!(first, second, "get_aligned_og must be deterministic");
    }

    #[test]
    fn test_get_aligned_og_large_window_does_not_panic() {
        let engine = simple_visualize_engine();
        // window larger than the genome — should clamp, not panic
        let result = engine.get_aligned_og(0, 1000, None);
        assert_eq!(result.len(), 2);
    }
}
