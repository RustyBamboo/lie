use approx::assert_abs_diff_eq;
use lie::gellmann::*;
use lie::lie_algebra::*;

#[test]
fn check_pauli_structure() {
    let g = gen_gellmann(2);
    let x = find_structure_constants(&g);

    // The Pauli matrices satisfy:
    // [sigma_x, sigma_y] = 2sigma_z

    assert_abs_diff_eq!(x.get(&(0usize, 1usize)).unwrap().1.re, 2., epsilon = 1e-8);
    assert_abs_diff_eq!(x.get(&(1usize, 2usize)).unwrap().1.re, 2., epsilon = 1e-8);
    assert_abs_diff_eq!(x.get(&(2usize, 0usize)).unwrap().1.re, 2., epsilon = 1e-8);
}

#[test]
fn test_gellmann_structure() {
    let g = gen_gellmann(3);
    let x = find_structure_constants(&g);

    println!("{}", x.len());
}
