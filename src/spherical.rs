use ndarray as nd;

use num_complex::Complex64;

use crate::su2;

pub fn c_g_p(r: i32, u: i32) -> f64 {
    let r = r as f64;
    let u = u as f64;
    let x: f64 = (((1. + r + u) * (2. + r + u)) / ((1. + r) * (1. + 2. * r))).sqrt();
    x / 2f64.sqrt()
}

pub fn c_g_z(r: i32, u: i32) -> f64 {
    let r = r as f64;
    let u = u as f64;
    (((1. + r - u) * (1. + r + u)) / ((1. + r) * (1. + 2. * r))).sqrt()
}

pub fn c_g_n(r: i32, u: i32) -> f64 {
    let r = r as f64;
    let u = u as f64;
    f64::sqrt(
        (2. + 3. * r + f64::powf(r, 2.) - 3. * u - 2. * r * u + f64::powf(u, 2.))
            / (2. + 6. * r + 4. * f64::powf(r, 2.)),
    )
}

pub fn c_g(r: i32, u: i32, la: i32) -> f64 {
    match la {
        1 => c_g_p(r, u),
        0 => c_g_z(r, u),
        -1 => c_g_n(r, u),
        _ => panic!("Input for Clebsh Gordan doesn't make sense!"),
    }
}

pub fn q_1_u(j: f64, u: i32) -> nd::Array2<f64> {
    use su2::{s_x, s_y, s_z};
    match u {
        0 => s_z(j),
        1 => -(s_x(j) + s_y(j)),
        -1 => (s_x(j) - s_y(j)),
        _ => panic!("Bad input for u"),
    }
}

pub fn q_r_u(j: f64, r: i32, u: i32) -> nd::Array2<f64> {
    if r == 1 {
        q_1_u(j, u)
    } else {
        let n = (j * 2. + 1.) as usize;

        let mut mat: nd::Array2<f64> = nd::Array2::zeros((n, n));

        for i in (-1 + u)..=(u + 1) {
            if i.abs() < r {
                let q1u = q_1_u(j, -i + u);
                mat = mat + c_g(r - 1, i, -i + u) * q_r_u(j, r - 1, i).dot(&q1u);
            }
        }
        mat
    }
}

pub fn basis_from_spin(j: f64) -> Vec<nd::Array2<f64>> {
    let n = (j * 2. + 1.) as u32;

    let mut basis = Vec::with_capacity(2usize.pow(n) - 1);

    for r in 1..(n as i32) {
        for u in -r..=r {
            let x = q_r_u(j, r, u);
            basis.push(x);
        }
    }
    basis
}

pub fn hermitian_basis_from_spin(j: f64) -> Vec<nd::Array2<Complex64>> {
    let n = (j * 2. + 1.) as u32;
    let basis = basis_from_spin(j);

    let mut herm_basis: Vec<nd::Array2<Complex64>> = Vec::with_capacity(2usize.pow(n) - 1);

    let mut i = 0;

    for r in 1..n {
        for u in 0..=r {
            let j: u32 = i + r + u;
            let mut sign = 1;
            if u % 2 == 1 {
                sign = -1;
            }
            let idx = j - 2 * u;
            let b = (sign as f64) * &basis[idx as usize];
            let u_a = &basis[j as usize];

            let u1 = (u_a + &b).map(|&x| Complex64::new(x, 0.));

            herm_basis.push(u1 / 2.);

            if u == 0 {
                continue;
            }

            let u2 = u_a - &b;
            herm_basis.push(u2.map(|&x| Complex64::new(0., x)) / 2.);
        }
        i += 2 * r + 1;
    }
    herm_basis
}
