# SymbolicFact.jl

This package is an optimized computing package for the symbolic factorization of sparse matrices. works both with square and rectangular matrices.

- compute the elimination tree of the matrix : [`etree`](src/etree.jl)
- compute the rows and columns count of the matrix, can detect simple supernodes : [`rowcolcount!`](src/rowcolcount!.jl)
- compute the row indices of the Cholesky factor (only for square matrices) : [`idx_r`](src/idx_r.jl)
- plot trees with few functions : [`treeplots.jl`](src/treeplots.jl)

To do all this as once, use the function [`symbolicfact`](src/SymbolicFact.jl).

## Example
We consider the following matrix:

$$
 M =
\begin{pmatrix} 
  X & . & X & X & . \\
  . & X & X & X & X \\
  X & X & X & . & X \\
  X & X & . & X & . \\
  . & X & X & . & X
\end{pmatrix}
$$

```julia
julia> M = sparse([1,1,1,2,2,2,2,3,3,3,3,4,4,4,5,5,5],[1,3,4,2,3,4,5,1,2,3,5,1,2,4,2,3,5],ones(Int,17))
5×5 SparseMatrixCSC{Int64, Int64} with 17 stored entries:
 1  ⋅  1  1  ⋅
 ⋅  1  1  1  1
 1  1  1  ⋅  1
 1  1  ⋅  1  ⋅
 ⋅  1  1  ⋅  1
```
Compute the elimination tree with [`etree`](src/etree.jl):

```julia
julia> tree = etree(M)
5-element Vector{Int64}:
 3
 3
 4
 5
 0
```
We can plot it : 
```julia
julia> treeplot(tree,path)
```
