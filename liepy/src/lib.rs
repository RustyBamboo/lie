use num_complex::Complex64;
use numpy::{PyArray2, PyReadonlyArray2, ToPyArray};
use pyo3::prelude::{pymodule, PyModule, PyResult, Python};

#[allow(unused_imports)]
use lie::lie_algebra::{
    cross, dot, find_d_coefficients, find_structure_constants, su_commutator,
};
#[allow(unused_imports)]
use lie::spherical::hermitian_basis_from_spin;

use lie::gellmann::get_gellmann;
use lie::sylvester::gen_sylvester;

#[pymodule]
fn su_genpy(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    #[pyfn(m, "gen_sylvester")]
    fn gen_sylvester_py<'py>(py: Python<'py>, d: usize) -> Vec<&'py PyArray2<Complex64>> {
        let basis = gen_sylvester(d);
        let basis: Vec<&PyArray2<Complex64>> = basis.iter().map(|&x| x.into_pyarray(py)).collect();

        return basis;
    }

    #[pyfn(m, "gen_gellmann")]
    fn gen_gellmann_py<'py>(py: Python<'py>, d: usize) -> Vec<&'py PyArray2<Complex64>> {
        let basis = get_gellmann(d);
        let basis: Vec<&PyArray2<Complex64>> = basis.iter().map(|x| x.to_pyarray(py)).collect();

        return basis;
    }

    #[pyfn(m, "get_structure_constants")]
    fn get_structure_constants_py<'py>(
        _py: Python<'py>,
        basis: Vec<PyReadonlyArray2<'py, Complex64>>,
    ) -> std::collections::HashMap<(usize, usize), (usize, Complex64)> {
        let basis = &basis.iter().map(|x| x.as_array().to_owned()).collect();

        let struct_consts = find_structure_constants(basis);

        return struct_consts;
    }

    #[pyfn(m, "get_d_coefficients")]
    fn get_d_coefficients_py<'py>(
        _py: Python<'py>,
        basis: Vec<PyReadonlyArray2<'py, Complex64>>,
    ) -> std::collections::HashMap<(usize, usize, usize), Complex64> {
        let basis = &basis.iter().map(|x| x.as_array().to_owned()).collect();

        let struct_consts = find_d_coefficients(basis);

        struct_consts
    }

    #[pyfn(m, "su_commutator")]
    fn su_commutator_py<'py>(
        py: Python<'py>,
        l_a: PyReadonlyArray2<'py, Complex64>,
        l_b: PyReadonlyArray2<'py, Complex64>,
        f_ijk: std::collections::HashMap<(usize, usize), (usize, Complex64)>,
        basis: Vec<PyReadonlyArray2<'py, Complex64>>,
    ) -> &'py PyArray2<Complex64> {
        let basis = &basis.iter().map(|x| x.as_array().to_owned()).collect();
        let res = su_commutator(
            &l_a.as_array().to_owned(),
            &l_b.as_array().to_owned(),
            &f_ijk,
            basis,
        );

        res.to_pyarray(py)
    }

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

    #[pyfn(m, "dot")]
    fn dot_py<'py>(
        py: Python<'py>,
        l_a: PyReadonlyArray2<'py, Complex64>,
        l_b: PyReadonlyArray2<'py, Complex64>,
        f_ijk: std::collections::HashMap<(usize, usize, usize), Complex64>,
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
