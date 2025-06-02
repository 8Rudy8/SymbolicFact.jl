# SymbolicFact.jl

This package is an optimized computing package for the symbolic factorization of sparse matrices. works both with square and rectangular matrices.

- compute the elimination tree of the matrix : [`etree`](src/etree.jl)
- compute the rows and columns count of the matrix, can detect simple supernodes : [`rowcolcount!`](src/rowcolcount!.jl)
- compute the row indices of the Cholesky factor (only for square matrices) : [`idx_r`](src/idx_r.jl)
- plot trees with few functions : [`treeplots.jl`](src/treeplots.jl)

To do all this as once, use the function [`symbolicfact`](src/SymbolicFact.jl).
