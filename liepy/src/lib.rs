use ndarray as nd;

use num_complex::Complex64;
use numpy::{PyArray2, PyReadonlyArray2, ToPyArray};
use pyo3::prelude::{pymodule, PyModule, PyResult, Python};

#[allow(unused_imports)]
use lie::lie_algebra::{
    cross, dot, find_d_coefficients, find_structure_constants, su_anticommutator, su_commutator,
};
#[allow(unused_imports)]
use lie::spherical::hermitian_basis_from_spin;

use lie::gellmann::gen_gellmann;
use lie::sylvester::gen_sylvester;

use lie::su2::gen_sl2;
use lie::su2::gen_su2;

#[pymodule]
fn liepy(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    /// Generate matrix representation of su(d) via clock and shift
    #[pyfn(m, "gen_sylvester")]
    fn gen_sylvester_py<'py>(py: Python<'py>, d: usize) -> Vec<&'py PyArray2<Complex64>> {
        let basis = gen_sylvester(d);
        let basis: Vec<&PyArray2<Complex64>> = basis.iter().map(|x| x.to_pyarray(py)).collect();

        basis
    }

    /// Generate matrix representation of su(d) via Generalized Gellmann method
    #[pyfn(m, "gen_gellmann")]
    fn gen_gellmann_py<'py>(py: Python<'py>, d: usize) -> Vec<&'py PyArray2<Complex64>> {
        let basis = gen_gellmann(d);
        let basis: Vec<&PyArray2<Complex64>> = basis.iter().map(|x| x.to_pyarray(py)).collect();

        basis
    }

    /// Generate matrix representation of su(2) for spin j
    #[pyfn(m, "gen_su2")]
    fn gen_su_py<'py>(py: Python<'py>, j: f64) -> Vec<&'py PyArray2<Complex64>> {
        let basis: Vec<&PyArray2<Complex64>> =
            gen_su2(j).iter().map(|x| x.to_pyarray(py)).collect();

        basis
    }

    /// Generate matrix representation of sl(2) for spin j
    #[pyfn(m, "gen_sl2")]
    fn gen_sl_py<'py>(py: Python<'py>, j: f64) -> Vec<&'py PyArray2<f64>> {
        let basis: Vec<&PyArray2<f64>> = gen_sl2(j).iter().map(|x| x.to_pyarray(py)).collect();

        basis
    }

    /// Find the commutation coefficients for matrices in su(d)
    #[pyfn(m, "get_structure_constants")]
    fn get_structure_constants_py<'py>(
        _py: Python<'py>,
        basis: Vec<PyReadonlyArray2<'py, Complex64>>,
    ) -> std::collections::HashMap<(usize, usize), (usize, Complex64)> {
        let basis: Vec<nd::Array2<Complex64>> =
            basis.iter().map(|x| x.as_array().to_owned()).collect();

        let struct_consts = find_structure_constants(basis.as_slice());

        struct_consts
    }

    /// Find the anti-commutation coefficients for matrices in su(d)
    #[pyfn(m, "get_d_coefficients")]
    fn get_d_coefficients_py<'py>(
        _py: Python<'py>,
        basis: Vec<PyReadonlyArray2<'py, Complex64>>,
    ) -> std::collections::HashMap<(usize, usize), (usize, Complex64)> {
        let basis: Vec<nd::Array2<Complex64>> =
            basis.iter().map(|x| x.as_array().to_owned()).collect();

        let struct_consts = find_d_coefficients(basis.as_slice());

        struct_consts
    }

    /// Find the commutation result of two matrices, given the structure constants
    #[pyfn(m, "su_commutator")]
    fn su_commutator_py<'py>(
        py: Python<'py>,
        l_a: PyReadonlyArray2<'py, Complex64>,
        l_b: PyReadonlyArray2<'py, Complex64>,
        f_ijk: std::collections::HashMap<(usize, usize), (usize, Complex64)>,
        basis: Vec<PyReadonlyArray2<'py, Complex64>>,
    ) -> &'py PyArray2<Complex64> {
        let basis: Vec<nd::Array2<Complex64>> =
            basis.iter().map(|x| x.as_array().to_owned()).collect();
        let res = su_commutator(
            &l_a.as_array().to_owned(),
            &l_b.as_array().to_owned(),
            &f_ijk,
            basis.as_slice(),
        );

        res.to_pyarray(py)
    }

    /// Find the anticommutation result of two matrices, given the structure constants
    #[pyfn(m, "su_anticommutator")]
    fn su_anticommutator_py<'py>(
        py: Python<'py>,
        l_a: PyReadonlyArray2<'py, Complex64>,
        l_b: PyReadonlyArray2<'py, Complex64>,
        d_ijk: std::collections::HashMap<(usize, usize), (usize, Complex64)>,
        basis: Vec<PyReadonlyArray2<'py, Complex64>>,
    ) -> &'py PyArray2<Complex64> {
        let basis: Vec<nd::Array2<Complex64>> =
            basis.iter().map(|x| x.as_array().to_owned()).collect();
        let res = su_anticommutator(
            &l_a.as_array().to_owned(),
            &l_b.as_array().to_owned(),
            &d_ijk,
            basis.as_slice(),
        );

        res.to_pyarray(py)
    }

    /// Compute the cross-product of two matrices of su(d), given the structure constants
    #[pyfn(m, "cross")]
    fn cross_py<'py>(
        py: Python<'py>,
        l_a: PyReadonlyArray2<'py, Complex64>,
        l_b: PyReadonlyArray2<'py, Complex64>,
        f_ijk: std::collections::HashMap<(usize, usize), (usize, Complex64)>,
    ) -> &'py PyArray2<Complex64> {
        let res = cross(
            &l_a.as_array().to_owned(),
            &l_b.as_array().to_owned(),
            &f_ijk,
        );

        res.to_pyarray(py)
    }

    /// Compute the dot-product of two matrices of su(d), given the structure constants
    #[pyfn(m, "dot")]
    fn dot_py<'py>(
        py: Python<'py>,
        l_a: PyReadonlyArray2<'py, Complex64>,
        l_b: PyReadonlyArray2<'py, Complex64>,
        f_ijk: std::collections::HashMap<(usize, usize), (usize, Complex64)>,
    ) -> &'py PyArray2<Complex64> {
        let res = dot(
            &l_a.as_array().to_owned(),
            &l_b.as_array().to_owned(),
            &f_ijk,
        );

        res.to_pyarray(py)
    }

    Ok(())
}
