![](https://img.shields.io/crates/v/lie?style=flat-square)
![](https://img.shields.io/crates/d/lie?label=crates.io%20downloads)
![](https://img.shields.io/pypi/pyversions/liepy?style=flat-square)
![](https://img.shields.io/pypi/wheel/liepy?style=flat-square)

# Lie

A numerical library for working with (representations) of Lie Groups and Algebras.

### Features

- [x] Spin matrices `su(2)` and ladder matrices `sl(2, C)`
- [x] Generalized [Gell-Mann Matrices](https://en.wikipedia.org/wiki/Gell-Mann_matrices) `su(n)` (Hermitian)
- [x] [Sylverster](https://en.wikipedia.org/wiki/Generalized_Clifford_algebra) "Clock" and "Shift" Matrices `su(n)` (non-Hermitian)
- [x] `su(n)` matrices by "Spherical Harmonics Tensors su(2) addition" via Clebsh-Gordan coefficients
- [x] Computation of [Structure constants](https://en.wikipedia.org/wiki/Structure_constants) 
- [x] Computation of dot/cross product in `su(n)` through structure constants

### Examples

#### Python bindings
```python
import liepy as lp
su_algebra = lp.gen_gellmann(3)
```
Additional examples, and notebooks, can found [here](liepy/examples).


#### Rust
```rust
use lie::gellmann::*;
use lie::lie_algebra::*;
let su_algebra = get_gellmann(3);
let f = find_structure_constants(su_algebra); 

println!("{:?}", f);
```

## Installation

Pre-built binary wheels are available. 

```
pip install liepy
```

## Compiling from source

`Lie` depends on:

- `openblas-devel`, e.g. (for Ubuntu/Debian `sudo apt install libopenblas-devel`) or equivalent
- [Rust](https://rustup.rs/) >= 1.58
- [Maturin](https://github.com/PyO3/maturin)

```
git clone https://github.com/RustyBamboo/lie
cd lie/liepy
maturin build --release --manylinux=off
pip3 install target/wheels/liepy-....whl --force-reinstall
```

### Building for manylinux

For manylinux compiled wheel, a Docker container is used.

```
cd lie
docker build -t maturin liepy/
docker run --rm -v $(pwd):/io maturin build --release -m liepy/Cargo.toml
```

## Tests

To ensure the library is working as intended, a test can be run:

```
cargo test
```

## License

`Lie` is free and opensource, released under MIT license.

