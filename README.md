# HOLZ: High-Order Entropy Encoding of Lempel-Ziv Factor Distances

Computes the High-Order Entropy Encoding for

 - the Lempel-Ziv 77 parsing
 - the bit-optimal parsing [1]

## build

```console
git clone https://github.com/nicolaprezza/high-order-lz77
git submodule init
git submodule update
mkdir build
cd build
cmake ..
make
```

## Thanks

We used [miltl](https://github.com/MitI-7)'s dynamic wavelet-matrix implementation for the experiments,
of which we include a modified version in the folder `WaveletMatrix`

## References
- Paolo Ferragina, Igor Nitto, Rossano Venturini: On the bit-complexity of Lempel-Ziv compression. SODA 2009: 768-777, http://dl.acm.org/citation.cfm?id=1496770.1496854

