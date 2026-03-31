use numpy::PyReadonlyArray1;
use pyo3::prelude::*;
use std::collections::HashMap;

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
    m.add_function(wrap_pyfunction!(calculate_synteny_ratio_py, m)?)?;
    Ok(())
}
