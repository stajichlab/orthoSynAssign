mod engine;
mod utils;

use crate::engine::*;
use crate::utils::*;
use pyo3::prelude::*;

#[pyfunction]
#[pyo3(name = "get_window")]
#[pyo3(text_signature = "(seqid_vec, og_vec, gene_idx, target_ogs, window_size, is_circular)")]
pub fn get_window_py(
    seqid_vec: Vec<i16>,
    og_vec: Vec<i32>,
    gene_idx: usize,
    mut target_ogs: Vec<i32>,
    window_size: usize,
    is_circular: bool,
) -> PyResult<Vec<usize>> {
    target_ogs.sort_unstable(); // Sort once for binary search
    let mut result_vec = Vec::with_capacity(window_size);

    get_window(
        &seqid_vec,
        &og_vec,
        gene_idx,
        &target_ogs,
        window_size,
        is_circular,
        &mut result_vec,
    );

    Ok(result_vec)
}

#[pyfunction]
#[pyo3(name = "calculate_synteny_ratio")]
#[pyo3(text_signature = "(win_a, win_b)")]
pub fn calculate_synteny_ratio_py(mut win_a: Vec<i32>, mut win_b: Vec<i32>) -> PyResult<f64> {
    // 1. Sort them here so we can use the slice-based logic
    win_a.sort_unstable();
    win_b.sort_unstable();

    // 2. Call the maintainable slice-based helper
    Ok(calculate_synteny_ratio(&win_a, &win_b))
}

#[pymodule]
fn rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<SyntenyEngine>()?;
    m.add_class::<VisualizeEngine>()?;
    m.add_function(wrap_pyfunction!(get_window_py, m)?)?;
    m.add_function(wrap_pyfunction!(calculate_synteny_ratio_py, m)?)?;
    Ok(())
}
