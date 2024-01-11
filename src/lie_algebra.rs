use ndarray as nd;

use ndarray_linalg::*;

use std::collections::HashMap;

use itertools::Itertools;

///
/// The structure constants of a lie algebra completely specify the product (commutator bracket) of
/// the algebra.
///
/// For example, any two element $T_a$ and $T_b$, the structure constant $f_{abc}$ defines the
/// relationship:
///
/// $$[T_a, T_b] = i \sum_{c=1}^{n^2 - 1} f_{abc} T_c$$
///
/// This function solves for $f_{abc}$ and returns a hashmap, with a key (a, b) and value (c,
/// $f_{abc}$. If $f_{abc} = 0$ then the HashMap will not have an entry.
///
pub fn find_structure_constants(
    basis: &[nd::Array2<c64>],
) -> HashMap<(usize, usize), (usize, c64)> {
    use approx::AbsDiffEq;
    use std::iter::FromIterator;

    let n_dim = basis[0].shape()[0];
    let n = n_dim * n_dim;

    let commutator = |x: &nd::Array2<c64>, y: &nd::Array2<c64>| x.dot(y) - y.dot(x);

    // Rewrite the basis as a matrix with each su(n) element mapped to a column
    // This then defines the matrix A in Ax=b
    let matrix: Vec<_> = basis
        .iter()
        .map(|x| nd::ArrayView::from_shape(n, x.as_slice().unwrap()).unwrap())
        .collect();
    let matrix = nd::stack(nd::Axis(1), matrix.as_slice()).unwrap();

    // HashMap to store the strucure constants
    let mut struct_consts = HashMap::new();

    for (i, t_a) in basis.iter().enumerate() {
        for (j, t_b) in basis.iter().enumerate() {
            // Ignore when same, since we know they commute
            if i == j {
                continue;
            }

            let f_t_c = commutator(t_a, t_b);

            let f_t_c = nd::Array::from_iter(f_t_c.iter().cloned());

            // TODO: This is dumb. Perhaps use QR decomposition to solve Ax=b where is A is not
            // square matrix. Maybe utilize the commutator relationship to simplify work
            let x = matrix.least_squares(&f_t_c).unwrap();
            let x = x.solution;

            let idx = x.iter().map(|x| x.abs_diff_ne(&c64::new(0., 0.), 1e-8));

            if idx.clone().filter(|x| *x).count() > 1usize {
                assert!(true, "This algebra does not satisfiey the lie algebra");
            }

            let idx = idx.into_iter().position(|x| x);

            if let Some(idx) = idx {
                struct_consts.insert((i, j), (idx, x[idx]));
            }
        }
    }

    struct_consts
}

pub fn find_d_coefficients(basis: &[nd::Array2<c64>]) -> HashMap<(usize, usize), (usize, c64)> {
    use approx::AbsDiffEq;
    use std::iter::FromIterator;

    let n_dim = basis[0].shape()[0];
    let n = n_dim * n_dim;

    let anti_commutator = |x: &nd::Array2<c64>, y: &nd::Array2<c64>| x.dot(y) + y.dot(x);

    let eye: nd::Array2<c64> = nd::Array2::eye(n_dim);

    // Rewrite the basis as a matrix with each su(n) element mapped to a column
    // This then defines the matrix A in Ax=b
    let matrix: Vec<_> = basis
        .iter()
        .map(|x| nd::ArrayView::from_shape(n, x.as_slice().unwrap()).unwrap())
        .collect();
    let matrix = nd::stack(nd::Axis(1), matrix.as_slice()).unwrap();

    // HashMap to store the strucure constants
    let mut struct_consts = HashMap::new();

    for (i, t_a) in basis.iter().enumerate() {
        for (j, t_b) in basis.iter().enumerate() {
            let mut f_t_c = anti_commutator(t_a, t_b);

            // if i == j {
            // f_t_c = f_t_c - &eye / c64::new(n_dim as f64, 0.);
            // }

            let f_t_c = nd::Array::from_iter(f_t_c.iter().cloned());

            // TODO: This is dumb. Perhaps use QR decomposition to solve Ax=b where is A is not
            // square matrix. Maybe utilize the commutator relationship to simplify work
            let x = matrix.least_squares(&f_t_c).unwrap();
            let x = x.solution;

            let idx = x.iter().map(|x| x.abs_diff_ne(&c64::new(0., 0.), 1e-8));

            let idx = idx.into_iter().position(|x| x);

            if let Some(idx) = idx {
                struct_consts.insert((i, j), (idx, x[idx]));
            }

            // for (c, val) in idx.clone().into_iter().enumerate() {
            //     if val {
            //         struct_consts.insert((i, j, c), x[c]);
            //     }
            // }
        }
    }

    struct_consts
}

///
/// Assume a vector space for the su lie algebra, with a vector being defined by the basis. Find
/// what the commutator is using the coordinates of two matrices defined on the basis.
/// l_a: coordinates of first matrix
/// l_b: coordinates of second matrix
/// f_ijk: structure constants for the lie algebra
/// basis: the vector basis for the lie algebra
///
/// returns the commutator result [l_a, l_b]
///
pub fn su_commutator(
    l_a: &ndarray::Array2<c64>,
    l_b: &ndarray::Array2<c64>,
    f_ijk: &HashMap<(usize, usize), (usize, c64)>,
    basis: &[nd::Array2<c64>],
) -> nd::Array2<c64> {
    let n_dim = basis[0].shape()[0];
    let mut res: nd::Array2<c64> = nd::Array2::zeros((n_dim, n_dim));
    for (i_i, l_i) in l_a.iter().enumerate() {
        for (i_j, l_j) in l_b.iter().enumerate() {
            let (i_k, f) = match f_ijk.get(&(i_i, i_j)) {
                Some(x) => x,
                None => continue,
            };
            let coeff = l_i * l_j * f;

            if coeff.abs() > 1e-8 {
                let prod = coeff * basis.get(*i_k).unwrap();
                res = res + prod;
            }
        }
    }
    res
}

///
/// Assume a vector space for the su lie algebra, with a vector being defined by the basis. Find
/// what the commutator is using the coordinates of two matrices defined on the basis.
/// l_a: coordinates of first matrix
/// l_b: coordinates of second matrix
/// f_ijk: structure constants for the lie algebra
/// basis: the vector basis for the lie algebra
///
/// returns the commutator result [l_a, l_b]
///
pub fn su_anticommutator(
    l_a: &ndarray::Array2<c64>,
    l_b: &ndarray::Array2<c64>,
    d_ijk: &HashMap<(usize, usize), (usize, c64)>,
    basis: &[nd::Array2<c64>],
) -> nd::Array2<c64> {
    let n_dim = basis[0].shape()[0];
    let mut res: nd::Array2<c64> = nd::Array2::zeros((n_dim, n_dim));
    for (i_i, l_i) in l_a.iter().enumerate() {
        for (i_j, l_j) in l_b.iter().enumerate() {
            let (i_k, f) = match d_ijk.get(&(i_i, i_j)) {
                Some(x) => x,
                None => continue,
            };
            let coeff = l_i * l_j * f;

            if coeff.abs() > 1e-8 {
                let prod = coeff * basis.get(*i_k).unwrap();
                res = res + prod;
            }
        }
    }
    res
}

pub fn cross(
    l_a: &ndarray::Array2<c64>,
    l_b: &ndarray::Array2<c64>,
    f_ijk: &HashMap<(usize, usize), (usize, c64)>,
) -> nd::Array2<c64> {
    let n_dim = l_a.shape()[0];
    let mut res: nd::Array2<c64> = nd::Array2::zeros((l_a.shape()[0], l_a.shape()[1]));
    for (i_i, l_i) in l_a.iter().enumerate() {
        for (i_j, l_j) in l_b.iter().enumerate() {
            for k in 0..n_dim {
                let (i_k, f) = match f_ijk.get(&(i_i, i_j)) {
                    Some(x) => x,
                    None => continue,
                };
                if k != *i_k {
                    continue;
                }

                res[[k, 0]] += f * l_i * l_j;
            }
        }
    }
    res
}

pub fn dot(
    l_a: &ndarray::Array2<c64>,
    l_b: &ndarray::Array2<c64>,
    d_ijk: &HashMap<(usize, usize), (usize, c64)>,
) -> nd::Array2<c64> {
    let n_dim = l_a.shape()[0];
    let mut res: nd::Array2<c64> = nd::Array2::zeros((l_a.shape()[0], l_a.shape()[1]));
    for (i_i, l_i) in l_a.iter().enumerate() {
        for (i_j, l_j) in l_b.iter().enumerate() {
            for k in 0..n_dim {
                let (i_k, f) = match d_ijk.get(&(i_i, i_j)) {
                    Some(x) => x,
                    None => continue,
                };
                if k != *i_k {
                    continue;
                }

                res[[k, 0]] += f * l_i * l_j;
            }
        }
    }
    res
}
