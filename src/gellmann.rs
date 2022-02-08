use ndarray as nd;
use num_complex::Complex64;

///
/// Construct an element for Generalized Gell-Mann Matrix
/// https://en.wikipedia.org/wiki/Generalizations_of_Pauli_matrices
///
pub fn gellmann(j: usize, k: usize, d: usize) -> nd::Array2<Complex64> {
    let mut gjkd: nd::Array2<Complex64> = nd::Array::zeros((d, d));

    let r1 = Complex64::new(1., 0.);
    let i1 = Complex64::new(0., 1.);

    if j > k {
        gjkd[[j - 1, k - 1]] = r1;
        gjkd[[k - 1, j - 1]] = r1;
    } else if k > j {
        gjkd[[j - 1, k - 1]] = -i1;
        gjkd[[k - 1, j - 1]] = i1;
    } else if j == k && j < d {
        let f = |n| {
            let n = n + 1;
            if n <= j {
                r1
            } else if n == (j + 1) {
                Complex64::new(-(j as f64), 0.)
            } else {
                Complex64::new(0., 0.)
            }
        };
        let diag = nd::Array1::from_shape_fn(d, f);
        let norm = (2. / (j as f64 * (j as f64 + 1.))).sqrt();
        let norm = Complex64::new(norm, 0.);

        gjkd = norm * nd::Array2::from_diag(&diag);
    } else {
        gjkd = nd::Array2::eye(d);
    }
    gjkd
}

///
/// Returns a basis for su(d) via Generalized Gell-Mann matrices
///
pub fn gen_gellmann(d: usize) -> Vec<nd::Array2<Complex64>> {
    let mut basis: Vec<nd::Array2<Complex64>> = Vec::with_capacity(d.pow(2));

    for j in 1..d + 1 {
        for k in 1..d + 1 {
            // Don't include identity
            if j == d && k == d {
                continue;
            }
            let el = Complex64::new(0., 1.) * gellmann(j, k, d);

            basis.push(el);
        }
    }
    basis
}
