<h1 align="center">arkworks::ldt</h1>

<p align="center">
    <a href="https://github.com/arkworks-rs/sponge/blob/master/LICENSE-APACHE">
        <img src="https://img.shields.io/badge/license-APACHE-blue.svg"></a>
    <a href="https://github.com/arkworks-rs/sponge/blob/master/LICENSE-MIT">
        <img src="https://img.shields.io/badge/license-MIT-blue.svg"></a>
</p>

`ark-ldt` is a Rust library that provides implementations of low-degree tests (LDTs). This library is released under the MIT License
and the Apache v2 License (see [License](#license)).

**WARNING:** This is an academic prototype, and in particular has not received careful code review.
This implementation is NOT ready for production use.

## Overview

A (univariate) low-degree test is an IOP that checks that a given function is close to a (univariate) polynomial of low degree. This library provides two LDTs: the **direct low-degree test** and the **FRI Protocol**. The library also comes with R1CS constraints for the LDT verifiers. Enable `r1cs` feature to use those constraints. 

## Build Guide

The library compiles on the `stable` toolchain of the Rust compiler. To install the latest version
of Rust, first install `rustup` by following the instructions [here](https://rustup.rs/), or via
your platform's package manager. Once `rustup` is installed, install the Rust toolchain by invoking:
```bash
rustup install stable
```

After that, use `cargo` (the standard Rust build tool) to build the library:
```bash
git clone https://github.com/arkworks-rs/ldt.git
cd sponge 
cargo build --release
```

This library comes with some unit and integration tests. Run these tests with:
```bash
cargo test
```


## License

This library is licensed under either of the following licenses, at your discretion.

* [Apache License Version 2.0](LICENSE-APACHE)
* [MIT License](LICENSE-MIT)

Unless you explicitly state otherwise, any contribution that you submit to this library shall be
dual licensed as above (as defined in the Apache v2 License), without any additional terms or
conditions.

## Reference papers

[Fractal: Post-Quantum and Transparent Recursive Proofs from Holography][cos20]<br>
Alessandro Chiesa, Dev Ojha, Nicholas Spooner     

[Fast Reed-Solomon Interactive Oracle Proofs of Proximity][bbhr17]<br>
Eli Ben-Sasson, Iddo Bentov, Ynon Horesh, Michael Riabzev

[cos20]: https://eprint.iacr.org/2019/1076
[bbhr17]: https://eccc.weizmann.ac.il/report/2017/134/
