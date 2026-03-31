use numpy::{IntoPyArray, PyArray1, PyReadonlyArray1};
use pyo3::prelude::*;
use std::collections::HashMap;

fn get_window_logic(
    og_masked: &[bool],
    seq: &[i16],
    gene_idx: usize,
    window_size: usize,
) -> Vec<usize> {
    let half_win = window_size / 2;
    let focal_seqid = seq[gene_idx];

    let valid_indices: Vec<usize> = (0..og_masked.len())
        .filter(|&i| og_masked[i] && seq[i] == focal_seqid)
        .collect();

    let pos = valid_indices.binary_search(&gene_idx).unwrap_or_else(|p| p);

    let left_start = pos.saturating_sub(half_win);
    let left_side = &valid_indices[left_start..pos];

    let start_r = if pos < valid_indices.len() && valid_indices[pos] == gene_idx {
        pos + 1
    } else {
        pos
    };

    let right_end = (start_r + half_win).min(valid_indices.len());
    let right_side = &valid_indices[start_r..right_end];

    // Combine them into a single Vec
    let mut result = Vec::with_capacity(left_side.len() + right_side.len());
    result.extend_from_slice(left_side);
    result.extend_from_slice(right_side);
    result
}

#[pyfunction]
#[pyo3(name = "get_window")]
pub fn get_window_py<'py>(
    py: Python<'py>,
    og_masked_array: PyReadonlyArray1<'py, bool>,
    seq_array: PyReadonlyArray1<'py, i16>,
    gene_idx: usize,
    window_size: usize,
) -> PyResult<Bound<'py, PyArray1<usize>>> {
    // Convert NumPy to Rust Slices
    let og_masked = og_masked_array.as_slice()?;
    let seq = seq_array.as_slice()?;

    // Call the pure Rust logic
    let result_vec = get_window_logic(og_masked, seq, gene_idx, window_size);

    // Convert Rust Vec back to NumPy
    // Since we aren't optimizing with 'unsafe' here for simplicity:
    Ok(result_vec.into_pyarray(py))
}

fn calculate_synteny_ratio_logic(win_a: &[i32], win_b: &[i32]) -> f64 {
    let len_a = win_a.len();
    let len_b = win_b.len();

    // Quick exit for empty windows
    if len_a == 0 || len_b == 0 {
        return 0.0;
    }

    // Count frequencies in win_a
    // For small windows, HashMap is great.
    // If windows were HUGE, we might sort and peek, but for n=20, this is ideal.
    let mut counts_a = HashMap::with_capacity(len_a);
    for &og in win_a {
        *counts_a.entry(og).or_insert(0) += 1;
    }

    // Count frequencies in win_b
    let mut counts_b = HashMap::with_capacity(len_b);
    for &og in win_b {
        *counts_b.entry(og).or_insert(0) += 1;
    }

    // Calculate matches (min count of shared IDs)
    let mut matches = 0;

    // Always iterate over the smaller map to minimize lookups
    let (small, large) = if counts_a.len() < counts_b.len() {
        (&counts_a, &counts_b)
    } else {
        (&counts_b, &counts_a)
    };

    for (og_id, count_small) in small {
        if let Some(count_large) = large.get(og_id) {
            matches += std::cmp::min(*count_small, *count_large);
        }
    }

    // Normalize by the longer window length
    matches as f64 / std::cmp::max(len_a, len_b) as f64
}

#[pyfunction]
#[pyo3(name = "calculate_synteny_ratio")]
pub fn calculate_synteny_ratio_py(
    win_a_array: PyReadonlyArray1<i32>,
    win_b_array: PyReadonlyArray1<i32>,
) -> PyResult<f64> {
    // Convert to slices (Zero-copy)
    let win_a = win_a_array.as_slice()?;
    let win_b = win_b_array.as_slice()?;

    // Call pure Rust logic
    Ok(calculate_synteny_ratio_logic(win_a, win_b))
}

#[pymodule]
fn rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // This registers your function so it appears as rs.get_window()
    m.add_function(wrap_pyfunction!(get_window_py, m)?)?;
    m.add_function(wrap_pyfunction!(calculate_synteny_ratio_py, m)?)?;
    Ok(())
}
