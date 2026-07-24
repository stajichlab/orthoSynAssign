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

pub fn get_window<F, G>(
    seqid_vec: &[i16],
    gene_idx: usize,
    window_size: usize,
    is_circular: bool,
    og_is_valid: F,
    og_value: G,
    buffer: &mut Vec<usize>,
) where
    F: Fn(usize) -> bool,
    G: Fn(usize) -> i32,
{
    if is_circular {
        get_window_circular(
            seqid_vec,
            gene_idx,
            window_size,
            og_is_valid,
            og_value,
            buffer,
        );
    } else {
        get_window_linear(
            seqid_vec,
            gene_idx,
            window_size,
            og_is_valid,
            og_value,
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

fn get_window_linear<F, G>(
    seqid_vec: &[i16],
    gene_idx: usize,
    window_size: usize,
    og_is_valid: F,
    og_value: G,
    buffer: &mut Vec<usize>,
) where
    F: Fn(usize) -> bool,
    G: Fn(usize) -> i32,
{
    buffer.clear();
    let half_win = window_size / 2;
    let focal_seqid = seqid_vec[gene_idx];

    // Look Left: scan backwards, skipping consecutive tandem repeats
    let mut left_indices = Vec::with_capacity(half_win);
    let (mut i, mut left_count) = (gene_idx, 0);
    let mut last_og = og_value(gene_idx); // initialize with focal gene's OG
    while i > 0 && left_count < half_win {
        i -= 1;
        if seqid_vec[i] == focal_seqid && og_is_valid(i) {
            let og = og_value(i);
            if og != last_og {
                left_indices.push(i);
                left_count += 1;
                last_og = og;
            }
        }
    }
    // Reverse to restore ascending order
    buffer.extend(left_indices.into_iter().rev());

    // Look Right: scan forwards, skipping consecutive tandem repeats
    let (mut j, mut right_count) = (gene_idx, 0);
    last_og = og_value(gene_idx); // reset to focal gene's OG
    while j < seqid_vec.len() - 1 && right_count < half_win {
        j += 1;
        if seqid_vec[j] == focal_seqid && og_is_valid(j) {
            let og = og_value(j);
            if og != last_og {
                buffer.push(j);
                right_count += 1;
                last_og = og;
            }
        }
    }
}

fn get_window_circular<F, G>(
    seqid_vec: &[i16],
    gene_idx: usize,
    window_size: usize,
    og_is_valid: F,
    og_value: G,
    buffer: &mut Vec<usize>,
) where
    F: Fn(usize) -> bool,
    G: Fn(usize) -> i32,
{
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
    let mut left_indices = Vec::with_capacity(half_win);
    let (mut i, mut left_count) = (gene_idx, 0);
    let mut last_og = og_value(gene_idx); // initialize with focal gene's OG
    let mut l_steps = 0;

    while left_count < half_win && l_steps < max_neighbors {
        i = (i + n_genes - 1) % n_genes;
        l_steps += 1;
        if i == gene_idx {
            break;
        }
        if seqid_vec[i] == focal_seqid && og_is_valid(i) {
            let og = og_value(i);
            if og != last_og {
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
    last_og = og_value(gene_idx); // reset to focal gene's OG
    let mut r_steps = 0;

    while right_count < r_limit && r_steps < max_neighbors - l_steps {
        j = (j + 1) % n_genes;
        r_steps += 1;
        if j == gene_idx {
            break;
        }

        if seqid_vec[j] == focal_seqid && og_is_valid(j) {
            let og = og_value(j);
            if og != last_og {
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
