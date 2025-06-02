# SymbolicFact.jl

This package is an optimized computing package for the symbolic factorization of sparse matrices. works both with square and rectangular matrices.

- compute the elimination tree of the matrix : [`etree`](@ref)
- compute the rows and columns count of the matrix, can detect simple supernodes : [`rowcol:count`](@ref)
- compute the row indices of the Cholesky factor (only for square matrices) : [`idx_r`](@ref)
- plot trees with few functions : [`treeplots.jl`](@ref)

To do all this as once, use the function [`symbolicfact`](@ref).
