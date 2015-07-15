# Sparse2015

This is a C implementation of several storage formats for sparse matrices, as well as the implementation of the sparse matrix vector product (SpMV) for each format. The following formats are implemented so far:

- COO: Coordinate list.
- CSR: Compressed Sparse Row.
- CSC: Compressed Sparse Column.
- OSKI: Compressed Sparse Row format, using Berkeley implementation for SpMV.
- SSS: Symmetric Sparse Skyline.
- TJDS: Transposed Jagged Diagonal Storage.

The provided Makefile can be used to compile on *nix or on Windows using MinGW.
