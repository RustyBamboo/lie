use ndarray as nd;

/// Helper for create raising operator elements
fn c_p(j: f64, m: &nd::Array1<f64>) -> nd::Array1<f64> {
    ((j - m) * (m + j + 1.)).mapv(f64::sqrt)
}

/// Helper for create lowering operator elements
fn c_m(j: f64, m: &nd::Array1<f64>) -> nd::Array1<f64> {
    ((m + j) * (-m + j + 1.)).mapv(f64::sqrt)
}

/// Create a raising matrix for a spin j
pub fn j_p(j: f64) -> nd::Array2<f64> {
    let n = (j * 2. + 1.) as usize;

    let cp = c_p(j, &nd::Array1::range(j - 1., -j - 1., -1.));

    let mut mat = nd::Array2::zeros((n, n));
    for i in 0..n - 1 {
        mat[[i, i + 1]] = cp[i];
    }
    mat
}

/// Create a lowering matrix for a spin j
pub fn j_m(j: f64) -> nd::Array2<f64> {
    let n = (j * 2. + 1.) as usize;

    let cm = c_m(j, &nd::Array1::range(j, -j, -1.));

    let mut mat = nd::Array2::zeros((n, n));
    for i in 1..n {
        mat[[i, i - 1]] = cm[i - 1];
    }
    mat
}

/// Create a s_z matrix for spin j
pub fn s_z(j: f64) -> nd::Array2<f64> {
    let n = (j * 2. + 1.) as usize;
    let mut mat = nd::Array2::zeros((n, n));

    let diag = nd::Array1::range(j, -j - 1., -1.);

    for i in 0..n {
        mat[[i, i]] = diag[i];
    }
    mat
}

/// Create a s_x matrix for spin j
pub fn s_x(j: f64) -> nd::Array2<f64> {
    (j_p(j) + j_m(j)) * 0.5
}

/// Create a s_y matrix for spin j
/// Note: this returns a real matrix to represent the complex matrix
pub fn s_y(j: f64) -> nd::Array2<f64> {
    (j_p(j) - j_m(j)) * 0.5
}

use num_complex::Complex64;
/// Returns matrix representations of s_z, s_x, s_y for spin j
/// For convenience, these matrices are returned as Complex
pub fn gen_su2(j: f64) -> [nd::Array2<Complex64>; 3] {
    let _i = Complex64::new(0., 1.);
    let _1 = Complex64::new(1., 0.);

    let s_z = s_z(j).map(|&x| Complex64::new(x, 0.));
    let s_x = s_x(j).map(|&x| Complex64::new(x, 0.));
    let s_y = s_y(j).map(|&x| Complex64::new(x, 0.));
    [s_z, s_x, s_y]
}

/// Returns matrix representations of diagonal s_z, raising and lowering operators s_+, s_-
pub fn gen_sl2(j: f64) -> [nd::Array2<f64>; 3] {
    [s_z(j), j_p(j), j_m(j)]
}
