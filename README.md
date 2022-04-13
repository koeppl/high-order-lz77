# HOLZ: High-Order Entropy Encoding of Lempel-Ziv Factor Distances

Computes the High-Order Entropy Encoding for

 - the Lempel-Ziv 77 parsing
 - the bit-optimal parsing [1]

## build

```console
git clone https://github.com/koeppl/high-order-lz77.git
cd high-order-lz77
git submodule init
git submodule update
mkdir build
cd build
cmake ..
make
./holz -i Makefile
```

## Caveats
The code assumes that there is no `\0` byte in the input. One way would
be to remove all `\0` bytes, but then the decompression would not return
the original file. Another way would be to use an escaper such as defined in
https://github.com/koeppl/stringology-rust
In the following example, we apply this escaper to use the byte 254 as an escape 
character, so '254' becomes '254 254', and '0' becomes '254 1'.

```console
	cargo build
	cargo run --bin escape -- -i <INFILE> -e 254 -f '0' -t '1' -r -o <OUTFILE>
```
(Maybe call `cargo update -p clap_derive --precise 3.0.0-beta.2` if `cargo build` fails.)
Appending the parameter '-r' reverts the escaping. 


## Thanks

We used [MitI_7](https://github.com/MitI-7)'s dynamic wavelet-matrix implementation for the experiments,
of which we include a modified version in the folder `WaveletMatrix`

## References
- Paolo Ferragina, Igor Nitto, Rossano Venturini: On the bit-complexity of Lempel-Ziv compression. SODA 2009: 768-777, http://dl.acm.org/citation.cfm?id=1496770.1496854

