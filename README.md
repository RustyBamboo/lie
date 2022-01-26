# Lie

A numerical library for working with Lie Groups and Algebras.

Currently, this library was built for my studies in quantum mechanics and computing. At the moment, this means the primary focus are U and SU groups.

- [x] Generalized Gell-Mann Matrices (SU)
- [x] Sylverster "Clock" and "Shift" Matrices
- [x] su(n) matrices by "spin su(2) addition" via Clebsh-Gordan coefficients
- [ ] other groups? SO(n)? Sp(n)? 


```rust
use lie::gellmann::*;
use lie::lie_algebra::*;
let su_algebra = get_gellmann(3);
let f = find_structure_constants(su_algebra); 

println!("{:?}", f);
```

## Python Bindings

A convenient way to work with this library is through the python bindings.

```python
import liepy as lp
su_algebra = lp.get_gellmann(3);
```

### Building for yourself

The quickest way to install the python packages is through the following:

> If changed Cargo.toml to use `openblas-system` (default) then install openblas (e.g. `sudo apt install libopenblas-dev`)

> If changed Cargo.toml to use `openblas-static` then install fortran compiler (e.g. `sudo apt install gfortran`)

And clone this repo (e.g. `git clone https://github.com/RustyBamboo/lie`)

```
cargo install maturin
cd lie/liepy
maturin build --release --manylinux=off
pip3 install target/wheels/liepy-....whl --force-reinstall
```

where the path to the wheel file is filled out appropriately to the version python you have.

### Building for manylinux

> WIP: manylinux and OpenBLAS not happy.

This is for building a generic linux wheel, such as to publish on PyPi:

```
cd lie
docker build -t maturin liepy/
docker run --rm -v $(pwd):/io maturin build --release -m liepy/Cargo.toml
```
