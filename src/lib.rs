pub mod gellmann;
pub mod lie_algebra;
pub mod spherical;
pub mod su2;
pub mod sylvester;
pub mod universal;

#[cfg(test)]
mod tests {
    use crate::lie_algebra::*;
    use crate::spherical::*;
    use approx::AbsDiffEq;
    use num_complex::Complex64;

    use std::f64::consts::PI;

    #[test]
    fn test_sylv() {
        use crate::sylvester::*;

        let _x = gen_sylvester(3);
    }

    #[test]
    fn check_clifford_algebra() {
        let d = 3.;

        let j = (d - 1.) / 2.;
        let basis = hermitian_basis_from_spin(j);

        let algebra = Complex64::new(0., 2. * PI / d).exp();
        println!("{}", algebra);

        for (i, a) in basis.iter().enumerate() {
            for (j, b) in basis.iter().enumerate() {
                if i == j {
                    continue;
                }
                let prod_ab = a.dot(b);
                let prod_ba = b.dot(a);

                let approx = prod_ab.abs_diff_eq(&prod_ba, 1e-5);
                println!("{} {}\n{}\n{}\n{}\n=======", i, j, prod_ab, prod_ba, approx);
            }
        }
    }

    ///
    /// The su(n) lie algebra of the SU(n) lie group contains a finite set of elements A_k which
    /// are represented by Hermitian matrices iA_k = H_k for Hermitian matrix H_k = H_k^\dagger
    ///
    #[test]
    fn test_hermitian() {
        let basis = hermitian_basis_from_spin(1.);

        for b in &basis {
            let conj = b.t().mapv(|x| x.conj());
            let approx = b.abs_diff_eq(&conj, 1e-8);
            assert!(approx, "The matrices are not Hermitian");
        }
    }

    #[test]
    fn test_commutation() {
        let basis = hermitian_basis_from_spin(0.5);
        let x = find_structure_constants(&basis);
        assert_eq!(x.len(), 3 * 2);

        let basis = hermitian_basis_from_spin(1.);
        let x = find_structure_constants(&basis);
        assert_eq!(x.len(), 25 * 2);
    }
}
