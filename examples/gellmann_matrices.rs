use lie::gellmann::*;

fn main() {
    let arg = std::env::args().nth(1).expect("Dimension argument?");

    let dim = arg.parse::<usize>().ok().expect("Dimension is not usize");

    let x = gen_gellmann(dim);

    for x_i in &x {
        println!("{}", x_i);
    }
}
