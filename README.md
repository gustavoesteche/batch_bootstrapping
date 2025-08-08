# Bootstrapping Framework Implementations

In this repository there are implementations of a framework that defines a faster asymptotic way to perform bootstrapping in Fully Homomorphic Encryption (FHE).
The project includes different implementations and supporting material in C++, SageMath, and LaTeX.

## Structure

* **`src_cpp/`** – C++ implementation using [NTL](https://libntl.org/), designed to verify and benchmark the proposed framework in a high-performance environment. The presented make file compiles all the
.cpp files under the tests folder. 
* **`src_sage/`** – SageMath implementation for prototyping, mathematical validation, and experimentation.
* **`src_latex/`** – LaTeX documents containing theoretical descriptions, analysis, and results related to the framework.

## Overview

Bootstrapping is a key operation in FHE schemes that refreshes ciphertexts to reduce noise and allow unlimited computations on encrypted data. The framework implemented here follows the approach of improving the **asymptotic complexity** of bootstrapping by optimizing polynomial arithmetic, modular reductions, and automorphisms. Its worth to mention that this implementation still needs optimizations to reach the 
the theoretical requirements, some operations can be improved: polynomial multiplication, automorphisms and polynomial decompositions are examples.

## Requirements

### For `src_cpp`:

* C++11
* [NTL](https://libntl.org/)
* GMP

### For `src_sage`:

* SageMath 9.0 or newer

### For `src_latex`:

* LaTeX distribution (e.g., TeX Live)
* `biblatex` and other standard math packages

## License

This project is released under the GNU General Public License (GPL).