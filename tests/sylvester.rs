use lie::lie_algebra::*;
use lie::sylvester::*;

#[test]
fn test_sylvester_structure() {
    let basis = gen_sylvester(3);
    let x = find_structure_constants(&basis);

    println!("{}", x.len());
}
