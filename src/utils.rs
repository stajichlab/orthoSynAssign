use std::collections::{HashMap, HashSet};

pub fn get_orthogroups_vec(num_orthogroups: usize, ogs: &[Vec<i32>]) -> Vec<Vec<(usize, usize)>> {
    let mut orthogroups: Vec<Vec<(usize, usize)>> = vec![Vec::new(); num_orthogroups];
    for (genome_idx, og_vec) in ogs.iter().enumerate() {
        for (gene_idx, &og_idx) in og_vec.iter().enumerate() {
            if og_idx >= 0 {
                let og_idx_usize = og_idx as usize;
                if og_idx_usize < num_orthogroups {
                    orthogroups[og_idx_usize].push((genome_idx, gene_idx));
                }
            }
        }
    }
    orthogroups
}

pub fn build_shared_matrix(ogs: &[Vec<i32>]) -> Vec<Vec<Vec<i32>>> {
    let num_genomes = ogs.len();
    let genome_sets: Vec<HashSet<i32>> = ogs
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

pub fn align_windows(
    og_windows: Vec<((usize, usize), Vec<usize>)>,
) -> Vec<((usize, usize), Vec<Option<usize>>)> {
    if og_windows.is_empty() {
        return Vec::new();
    }

    let mut max_prefix_len = 0;
    let mut focal_gene_offsets = Vec::with_capacity(og_windows.len());

    // Get the indices of each focal gene in the Vec
    for (focal_gene, genes) in og_windows.iter() {
        let pos = genes.iter().position(|&x| x == focal_gene.1).unwrap_or(0);
        focal_gene_offsets.push(pos);
        if pos > max_prefix_len {
            max_prefix_len = pos;
        }
    }

    // Get the window with max total len after aligning the focal genes
    let mut max_total_len = 0;
    for (i, (_, genes)) in og_windows.iter().enumerate() {
        let total_len = (max_prefix_len - focal_gene_offsets[i]) + genes.len();
        if total_len > max_total_len {
            max_total_len = total_len;
        }
    }

    // Padding
    let mut aligned = Vec::with_capacity(og_windows.len());
    for ((focal_gene, genes), offset) in og_windows.into_iter().zip(focal_gene_offsets) {
        let front_pad_size = max_prefix_len - offset;
        let mut aligned_row = Vec::with_capacity(max_total_len);
        for _ in 0..front_pad_size {
            aligned_row.push(None);
        }
        aligned_row.extend(genes.into_iter().map(Some));

        while aligned_row.len() < max_total_len {
            aligned_row.push(None);
        }
        aligned.push((focal_gene, aligned_row));
    }
    aligned
}

pub fn get_window(
    seqid_vec: &[i16],
    ogs_vec: &[i32],
    gene_idx: usize,
    shared_ogs: &[i32],
    window_size: usize,
    is_circular: bool,
    buffer: &mut Vec<usize>,
) {
    if is_circular {
        get_window_circular(
            seqid_vec,
            ogs_vec,
            gene_idx,
            shared_ogs,
            window_size,
            buffer,
        );
    } else {
        get_window_linear(
            seqid_vec,
            ogs_vec,
            gene_idx,
            shared_ogs,
            window_size,
            buffer,
        );
    }
}

pub fn calculate_synteny_ratio(win_a: &[i32], win_b: &[i32]) -> f64 {
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

pub fn cluster_genes(
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
            let root_u = find_dsu(&mut dsu, u_id);
            let root_v = find_dsu(&mut dsu, v_id);

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
        let r = find_dsu(&mut dsu, i);
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

fn get_window_linear(
    seqid_vec: &[i16],
    ogs_vec: &[i32],
    gene_idx: usize,
    shared_ogs: &[i32],
    window_size: usize,
    buffer: &mut Vec<usize>,
) {
    buffer.clear();
    let half_win = window_size / 2;
    let focal_seqid = seqid_vec[gene_idx];

    // Look Left: scan backwards, skipping consecutive tandem repeats
    let mut last_og = ogs_vec[gene_idx]; // initialize with focal gene's OG
    let mut left_indices = Vec::with_capacity(half_win);
    let (mut i, mut left_count) = (gene_idx, 0);
    while i > 0 && left_count < half_win {
        i -= 1;
        if seqid_vec[i] == focal_seqid {
            let og = ogs_vec[i];
            if shared_ogs.binary_search(&og).is_ok() && og != last_og {
                left_indices.push(i);
                left_count += 1;
                last_og = og;
            }
        }
    }
    // Reverse to restore ascending order
    buffer.extend(left_indices.into_iter().rev());

    // Look Right: scan forwards, skipping consecutive tandem repeats
    last_og = ogs_vec[gene_idx]; // reset to focal gene's OG
    let (mut j, mut right_count) = (gene_idx, 0);
    while j < seqid_vec.len() - 1 && right_count < half_win {
        j += 1;
        if seqid_vec[j] == focal_seqid {
            let og = ogs_vec[j];
            if shared_ogs.binary_search(&og).is_ok() && og != last_og {
                buffer.push(j);
                right_count += 1;
                last_og = og;
            }
        }
    }
}

fn get_window_circular(
    seqid_vec: &[i16],
    ogs_vec: &[i32],
    gene_idx: usize,
    shared_ogs: &[i32],
    window_size: usize,
    buffer: &mut Vec<usize>,
) {
    buffer.clear();
    let half_win = window_size / 2;
    let focal_seqid = seqid_vec[gene_idx];
    let n_genes = seqid_vec.len();
    if n_genes == 0 {
        return;
    }
    let max_neighbors = std::cmp::min(window_size, n_genes - 1);

    // Look Left (Counter-clockwise)
    // l_steps always increments — hard cap against highly tandem-repeated short genomes
    let mut last_og = ogs_vec[gene_idx]; // initialize with focal gene's OG
    let mut left_indices = Vec::with_capacity(half_win);
    let (mut i, mut left_count) = (gene_idx, 0);
    let mut l_steps = 0;

    while left_count < half_win && l_steps < max_neighbors {
        i = (i + n_genes - 1) % n_genes;
        l_steps += 1;
        if i == gene_idx {
            break;
        }
        if seqid_vec[i] == focal_seqid {
            let og = ogs_vec[i];
            if shared_ogs.binary_search(&og).is_ok() && og != last_og {
                left_indices.push(i);
                left_count += 1;
                last_og = og;
            }
        }
    }
    buffer.extend(left_indices.into_iter().rev());

    // Look Right (Clockwise)
    // r_steps budget is what remains after left consumed l_steps
    let (mut j, mut right_count) = (gene_idx, 0);
    let r_limit = std::cmp::min(half_win, max_neighbors - l_steps);
    last_og = ogs_vec[gene_idx]; // reset to focal gene's OG
    let mut r_steps = 0;

    while right_count < r_limit && r_steps < max_neighbors - l_steps {
        j = (j + 1) % n_genes;
        r_steps += 1;
        if j == gene_idx {
            break;
        }

        if seqid_vec[j] == focal_seqid {
            let og = ogs_vec[j];
            if shared_ogs.binary_search(&og).is_ok() && og != last_og {
                buffer.push(j);
                right_count += 1;
                last_og = og;
            }
        }
    }
}

fn find_dsu(dsu: &mut [i32], mut i: usize) -> usize {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_orthogroups_vec_basic() {
        // 2 genomes, 3 orthogroups
        // genome 0: genes assigned to OGs [0, 1, 2]
        // genome 1: genes assigned to OGs [0, 2, -1]
        let ogs = vec![vec![0i32, 1, 2], vec![0, 2, -1]];
        let result = get_orthogroups_vec(3, &ogs);

        assert_eq!(result.len(), 3);
        assert_eq!(result[0], vec![(0, 0), (1, 0)]);
        assert_eq!(result[1], vec![(0, 1)]);
        assert_eq!(result[2], vec![(0, 2), (1, 1)]);
    }

    #[test]
    fn test_get_orthogroups_vec_negative_ignored() {
        let ogs = vec![vec![-1i32, -1, -1]];
        let result = get_orthogroups_vec(3, &ogs);
        assert_eq!(result.len(), 3);
        for og in &result {
            assert!(og.is_empty());
        }
    }

    #[test]
    fn test_get_orthogroups_vec_out_of_range_ignored() {
        // OG index 5 is >= num_orthogroups (3), should be ignored
        let ogs = vec![vec![0i32, 5, 2]];
        let result = get_orthogroups_vec(3, &ogs);
        assert_eq!(result[0], vec![(0, 0)]);
        assert!(result[1].is_empty());
        assert_eq!(result[2], vec![(0, 2)]);
    }

    #[test]
    fn test_get_orthogroups_vec_empty_ogs() {
        let ogs: Vec<Vec<i32>> = vec![];
        let result = get_orthogroups_vec(3, &ogs);
        assert_eq!(result.len(), 3);
        for og in &result {
            assert!(og.is_empty());
        }
    }

    #[test]
    fn test_get_orthogroups_vec_zero_orthogroups() {
        let ogs = vec![vec![0i32, 1]];
        let result = get_orthogroups_vec(0, &ogs);
        assert!(result.is_empty());
    }

    #[test]
    fn test_get_orthogroups_vec_empty_genome_vecs() {
        let ogs = vec![vec![], vec![]];
        let result = get_orthogroups_vec(2, &ogs);
        assert_eq!(result.len(), 2);
        for og in &result {
            assert!(og.is_empty());
        }
    }

    #[test]
    fn test_get_orthogroups_vec_multiple_genes_same_og_same_genome() {
        // Genome 0 has two genes both in OG 0
        let ogs = vec![vec![0i32, 0]];
        let result = get_orthogroups_vec(1, &ogs);
        assert_eq!(result[0], vec![(0, 0), (0, 1)]);
    }

    // --- build_shared_matrix ---

    #[test]
    fn test_build_shared_matrix_single_genome() {
        let ogs = vec![vec![1i32, 2, 3]];
        let matrix = build_shared_matrix(&ogs);
        assert_eq!(matrix.len(), 1);
        assert_eq!(matrix[0][0], vec![1, 2, 3]);
    }

    #[test]
    fn test_build_shared_matrix_two_genomes_full_overlap() {
        let ogs = vec![vec![1i32, 2, 3], vec![1, 2, 3]];
        let matrix = build_shared_matrix(&ogs);
        assert_eq!(matrix[0][1], vec![1, 2, 3]);
        assert_eq!(matrix[1][0], vec![1, 2, 3]);
    }

    #[test]
    fn test_build_shared_matrix_two_genomes_no_overlap() {
        let ogs = vec![vec![1i32, 2], vec![3i32, 4]];
        let matrix = build_shared_matrix(&ogs);
        assert!(matrix[0][1].is_empty());
        assert!(matrix[1][0].is_empty());
    }

    #[test]
    fn test_build_shared_matrix_two_genomes_partial_overlap() {
        let ogs = vec![vec![1i32, 2, 3], vec![2i32, 3, 4]];
        let matrix = build_shared_matrix(&ogs);
        assert_eq!(matrix[0][1], vec![2, 3]);
        assert_eq!(matrix[1][0], vec![2, 3]);
    }

    #[test]
    fn test_build_shared_matrix_ignores_minus_one() {
        let ogs = vec![vec![-1i32, 1, 2], vec![-1i32, 1, 3]];
        let matrix = build_shared_matrix(&ogs);
        // -1 should not appear in shared
        assert!(!matrix[0][1].contains(&-1));
        assert_eq!(matrix[0][1], vec![1]);
    }

    #[test]
    fn test_build_shared_matrix_diagonal_is_self() {
        let ogs = vec![vec![1i32, 2, 3], vec![4i32, 5]];
        let matrix = build_shared_matrix(&ogs);
        assert_eq!(matrix[0][0], vec![1, 2, 3]);
        assert_eq!(matrix[1][1], vec![4, 5]);
    }

    #[test]
    fn test_build_shared_matrix_symmetric() {
        let ogs = vec![vec![1i32, 2], vec![2i32, 3], vec![1i32, 3]];
        let matrix = build_shared_matrix(&ogs);
        for i in 0..3 {
            for j in 0..3 {
                assert_eq!(matrix[i][j], matrix[j][i]);
            }
        }
    }

    #[test]
    fn test_build_shared_matrix_sorted_output() {
        // Ensure the shared OGs are sorted
        let ogs = vec![vec![3i32, 1, 2], vec![2i32, 3, 1]];
        let matrix = build_shared_matrix(&ogs);
        let shared = &matrix[0][1];
        let mut sorted = shared.clone();
        sorted.sort_unstable();
        assert_eq!(shared, &sorted);
    }

    #[test]
    fn test_build_shared_matrix_empty_genomes() {
        let ogs = vec![vec![], vec![]];
        let matrix = build_shared_matrix(&ogs);
        assert!(matrix[0][1].is_empty());
    }

    // --- align_windows ---

    #[test]
    fn test_align_windows_empty() {
        let result = align_windows(vec![]);
        assert!(result.is_empty());
    }

    #[test]
    fn test_align_windows_single_row_focal_at_start() {
        // focal gene index is 0, gene list is [0, 1, 2]
        let input = vec![((0usize, 0usize), vec![0usize, 1, 2])];
        let result = align_windows(input);
        assert_eq!(result.len(), 1);
        let (focal, row) = &result[0];
        assert_eq!(*focal, (0, 0));
        assert_eq!(row, &vec![Some(0), Some(1), Some(2)]);
    }

    #[test]
    fn test_align_windows_single_row_focal_in_middle() {
        // focal gene index is 1, gene list is [0, 1, 2]
        let input = vec![((0usize, 1usize), vec![0usize, 1, 2])];
        let result = align_windows(input);
        assert_eq!(result.len(), 1);
        let (_, row) = &result[0];
        assert_eq!(row, &vec![Some(0), Some(1), Some(2)]);
    }

    #[test]
    fn test_align_windows_two_rows_same_offset() {
        // Both have focal gene at position 1
        let input = vec![
            ((0usize, 1usize), vec![0usize, 1, 2]),
            ((1usize, 4usize), vec![3usize, 4, 5]),
        ];
        let result = align_windows(input);
        assert_eq!(result.len(), 2);
        // No front padding needed since both offsets are equal
        assert_eq!(result[0].1, vec![Some(0), Some(1), Some(2)]);
        assert_eq!(result[1].1, vec![Some(3), Some(4), Some(5)]);
    }

    #[test]
    fn test_align_windows_front_padding() {
        // Row 0: focal at position 0 in [10, 11, 12]
        // Row 1: focal at position 2 in [20, 21, 22, 23]
        // max_prefix_len = 2
        // Row 0 gets 2 front Nones; Row 1 gets 0 front Nones
        let input = vec![
            ((0usize, 10usize), vec![10usize, 11, 12]),
            ((1usize, 22usize), vec![20usize, 21, 22, 23]),
        ];
        let result = align_windows(input);
        assert_eq!(result.len(), 2);
        let (_, row0) = &result[0];
        let (_, row1) = &result[1];
        assert_eq!(row0[0], None);
        assert_eq!(row0[1], None);
        assert_eq!(row0[2], Some(10));
        assert_eq!(row1[0], Some(20));
        assert_eq!(row1[2], Some(22));
    }

    #[test]
    fn test_align_windows_trailing_padding() {
        // Row 0: focal at position 1 in [0, 1, 2, 3, 4] (5 elements)
        // Row 1: focal at position 1 in [5, 6, 7] (3 elements)
        // Both have same prefix length; Row 1 gets trailing Nones
        let input = vec![
            ((0usize, 1usize), vec![0usize, 1, 2, 3, 4]),
            ((1usize, 6usize), vec![5usize, 6, 7]),
        ];
        let result = align_windows(input);
        let (_, row1) = &result[1];
        assert_eq!(row1.len(), 5);
        assert_eq!(row1[3], None);
        assert_eq!(row1[4], None);
    }

    #[test]
    fn test_align_windows_all_rows_same_length() {
        let input = vec![
            ((0usize, 0usize), vec![0usize, 1, 2]),
            ((1usize, 3usize), vec![3usize, 4, 5]),
            ((2usize, 6usize), vec![6usize, 7, 8]),
        ];
        let result = align_windows(input);
        let len = result[0].1.len();
        for (_, row) in &result {
            assert_eq!(row.len(), len);
        }
    }

    #[test]
    fn test_align_windows_focal_not_found_defaults_to_zero() {
        // focal gene index 99 not in genes list; unwrap_or(0) kicks in
        let input = vec![((0usize, 99usize), vec![0usize, 1, 2])];
        let result = align_windows(input);
        // Should not panic; offset defaults to 0
        assert_eq!(result.len(), 1);
    }

    #[test]
    fn test_align_windows_preserves_focal_gene_key() {
        let input = vec![((3usize, 7usize), vec![6usize, 7, 8])];
        let result = align_windows(input);
        assert_eq!(result[0].0, (3, 7));
    }

    #[test]
    fn test_get_window_linear_basic() {
        // Genome layout:
        // Index:      0   1   2   3   4
        // seqid_vec:  1,  1,  1,  1,  1
        // ogs_vec:   10, 20, 30, 40, 50
        let seqids = vec![1, 1, 1, 1, 1];
        let ogs = vec![10, 20, 30, 40, 50];
        let shared = vec![10, 20, 30, 40, 50]; // sorted for binary_search
        let mut buffer = Vec::new();

        // Focal gene at index 2 (og: 30), window size 4 (half_win = 2 left, 2 right)
        get_window(&seqids, &ogs, 2, &shared, 4, false, &mut buffer);

        // Expected left: [0, 1], right: [3, 4] -> combined: [0, 1, 3, 4]
        assert_eq!(buffer, vec![0, 1, 3, 4]);
    }

    #[test]
    fn test_calculate_synteny_ratio_empty_a() {
        assert_eq!(calculate_synteny_ratio(&[], &[1, 2, 3]), 0.0);
    }

    #[test]
    fn test_calculate_synteny_ratio_empty_b() {
        assert_eq!(calculate_synteny_ratio(&[1, 2, 3], &[]), 0.0);
    }

    #[test]
    fn test_calculate_synteny_ratio_both_empty() {
        assert_eq!(calculate_synteny_ratio(&[], &[]), 0.0);
    }

    #[test]
    fn test_calculate_synteny_ratio_identical() {
        let win = vec![1, 2, 3, 4, 5];
        assert_eq!(calculate_synteny_ratio(&win, &win), 1.0);
    }

    #[test]
    fn test_calculate_synteny_ratio_no_overlap() {
        let win_a = vec![1, 2, 3];
        let win_b = vec![4, 5, 6];
        assert_eq!(calculate_synteny_ratio(&win_a, &win_b), 0.0);
    }

    #[test]
    fn test_calculate_synteny_ratio_partial_overlap() {
        // win_a = [1, 2, 3], win_b = [2, 3, 4, 5]
        // matches = 2, max_len = 4
        let win_a = vec![1, 2, 3];
        let win_b = vec![2, 3, 4, 5];
        let ratio = calculate_synteny_ratio(&win_a, &win_b);
        assert!((ratio - 0.5).abs() < f64::EPSILON);
    }

    #[test]
    fn test_calculate_synteny_ratio_single_element_match() {
        let win_a = vec![42];
        let win_b = vec![42];
        assert_eq!(calculate_synteny_ratio(&win_a, &win_b), 1.0);
    }

    #[test]
    fn test_calculate_synteny_ratio_single_element_no_match() {
        let win_a = vec![1];
        let win_b = vec![2];
        assert_eq!(calculate_synteny_ratio(&win_a, &win_b), 0.0);
    }

    #[test]
    fn test_calculate_synteny_ratio_uses_max_len_as_denominator() {
        // win_a has 1 element, win_b has 10 elements, 1 match
        // ratio = 1/10
        let win_a = vec![5];
        let win_b = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        let ratio = calculate_synteny_ratio(&win_a, &win_b);
        assert!((ratio - 0.1).abs() < f64::EPSILON);
    }

    #[test]
    fn test_cluster_genes_empty_pairs() {
        let all_genes = vec![(0usize, 0usize), (0, 1), (1, 0)];
        let pairs = vec![];
        let result = cluster_genes(pairs, &all_genes);
        // Each gene is its own cluster
        assert_eq!(result.len(), 3);
        for cluster in &result {
            assert_eq!(cluster.len(), 1);
        }
    }

    #[test]
    fn test_cluster_genes_all_connected() {
        let all_genes = vec![(0usize, 0usize), (0, 1), (0, 2)];
        let pairs = vec![
            ((0usize, 0usize), (0usize, 1usize)),
            ((0usize, 1usize), (0usize, 2usize)),
        ];
        let result = cluster_genes(pairs, &all_genes);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].len(), 3);
    }

    #[test]
    fn test_cluster_genes_two_separate_clusters() {
        let all_genes = vec![(0, 0), (0, 1), (1, 0), (1, 1)];
        let pairs = vec![
            ((0usize, 0usize), (0usize, 1usize)),
            ((1usize, 0usize), (1usize, 1usize)),
        ];
        let result = cluster_genes(pairs, &all_genes);
        assert_eq!(result.len(), 2);
        assert_eq!(result[0].len(), 2);
        assert_eq!(result[1].len(), 2);
    }

    #[test]
    fn test_cluster_genes_self_pair_ignored() {
        let all_genes = vec![(0usize, 0usize), (0, 1)];
        let pairs = vec![((0usize, 0usize), (0usize, 0usize))];
        let result = cluster_genes(pairs, &all_genes);
        // Self-pair doesn't merge anything; 2 separate clusters
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_cluster_genes_unknown_gene_in_pair_ignored() {
        let all_genes = vec![(0usize, 0usize), (0, 1)];
        // (99, 99) not in all_genes
        let pairs = vec![((0usize, 0usize), (99usize, 99usize))];
        let result = cluster_genes(pairs, &all_genes);
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_cluster_genes_result_sorted() {
        let all_genes = vec![(1, 0), (0, 0), (0, 1)];
        let pairs = vec![((0usize, 0usize), (0usize, 1usize))];
        let result = cluster_genes(pairs, &all_genes);
        // Clusters sorted by first element; within cluster sorted too
        assert_eq!(result[0][0], (0, 0));
        assert_eq!(result[1][0], (1, 0));
    }

    #[test]
    fn test_cluster_genes_empty_genes() {
        let all_genes: Vec<(usize, usize)> = vec![];
        let pairs = vec![];
        let result = cluster_genes(pairs, &all_genes);
        assert!(result.is_empty());
    }

    #[test]
    fn test_get_window_linear_filters_tandem_and_unshared() {
        // Index:      0   1   2   3   4   5
        // seqid:      1,  1,  1,  1,  1,  1
        // ogs:       10, 10, 99, 20, 30, 40
        // - Index 0 & 1 are tandem repeats (og 10)
        // - Index 2 (og 99) is NOT in shared_ogs
        let seqids = vec![1, 1, 1, 1, 1, 1];
        let ogs = vec![10, 10, 99, 20, 30, 40];
        let shared = vec![10, 20, 30, 40]; // 99 missing
        let mut buffer = Vec::new();

        // Focal at index 3 (og 20), look left
        get_window(&seqids, &ogs, 3, &shared, 4, false, &mut buffer);

        // Left scan sees:
        // - Index 2: 99 (ignored, not shared)
        // - Index 1: 10 (kept)
        // - Index 0: 10 (ignored, tandem repeat of index 1)
        // Right scan sees:
        // - Index 4: 30 (kept)
        // - Index 5: 40 (kept)
        assert_eq!(buffer, vec![1, 4, 5]);
    }

    #[test]
    fn test_get_window_circular_wrap_around() {
        // Index:      0   1   2   3
        // seqids:     1,  1,  1,  1
        // ogs:       10, 20, 30, 40
        let seqids = vec![1, 1, 1, 1];
        let ogs = vec![10, 20, 30, 40];
        let shared = vec![10, 20, 30, 40];
        let mut buffer = Vec::new();

        // Focal at index 0 (og 10), window size 2 (1 left, 1 right)
        get_window(&seqids, &ogs, 0, &shared, 2, true, &mut buffer);

        // Left wraps to end (index 3), Right goes to index 1
        assert_eq!(buffer, vec![3, 1]);
    }

    #[test]
    fn test_get_window_ignores_different_seqids() {
        // Index:      0   1   2   3   4
        // seqids:     1,  2,  1,  2,  1  <-- mixed contigs/chromosomes
        // ogs:       10, 20, 30, 40, 50
        let seqids = vec![1, 2, 1, 2, 1];
        let ogs = vec![10, 20, 30, 40, 50];
        let shared = vec![10, 20, 30, 40, 50];
        let mut buffer = Vec::new();

        // Focal at index 2 (seqid = 1)
        get_window(&seqids, &ogs, 2, &shared, 4, false, &mut buffer);

        // Should only match indices with seqid == 1 (0 and 4)
        assert_eq!(buffer, vec![0, 4]);
    }

    #[test]
    fn test_find_dsu_single_element() {
        let mut dsu = vec![-1];
        assert_eq!(find_dsu(&mut dsu, 0), 0);
    }

    #[test]
    fn test_find_dsu_root_is_self() {
        // All elements are roots (negative values)
        let mut dsu = vec![-1, -1, -1];
        assert_eq!(find_dsu(&mut dsu, 0), 0);
        assert_eq!(find_dsu(&mut dsu, 1), 1);
        assert_eq!(find_dsu(&mut dsu, 2), 2);
    }

    #[test]
    fn test_find_dsu_chain() {
        // 0 -> 1 -> 2 (root)
        // dsu[0] = 1, dsu[1] = 2, dsu[2] = -1
        let mut dsu = vec![1, 2, -1];
        assert_eq!(find_dsu(&mut dsu, 0), 2);
        // Path compression: dsu[0] should now point directly to root (2)
        assert_eq!(dsu[0], 2);
    }

    #[test]
    fn test_find_dsu_path_compression() {
        // Chain: 0 -> 1 -> 2 -> 3 (root)
        let mut dsu = vec![1i32, 2, 3, -1];
        assert_eq!(find_dsu(&mut dsu, 0), 3);
        // After path compression, all nodes should point directly to root
        assert_eq!(dsu[0], 3);
        assert_eq!(dsu[1], 3);
        assert_eq!(dsu[2], 3);
        assert_eq!(dsu[3], -1); // root unchanged
    }

    #[test]
    fn test_find_dsu_already_compressed() {
        // 0 -> 2 (root), 1 -> 2 (root)
        let mut dsu = vec![2i32, 2, -1];
        assert_eq!(find_dsu(&mut dsu, 0), 2);
        assert_eq!(find_dsu(&mut dsu, 1), 2);
        assert_eq!(find_dsu(&mut dsu, 2), 2);
    }

    #[test]
    fn test_find_dsu_two_separate_trees() {
        // Tree 1: 0 -> 1 (root), Tree 2: 2 -> 3 (root)
        let mut dsu = vec![1i32, -1, 3, -1];
        assert_eq!(find_dsu(&mut dsu, 0), 1);
        assert_eq!(find_dsu(&mut dsu, 2), 3);
    }
}
