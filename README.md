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
julia> treeplot(tree,a_path)
347×179 Array{RGBA{N0f8},2} with eltype ColorTypes.RGBA{FixedPointNumbers.N0f8}:
 RGBA{N0f8}(1.0,1.0,1.0,1.0) ...
 ...
```
 ![elimination tree of M](/img/pippo.png)

Copmpute it's postorder and also plot it :
```julia
julia> porder = postorder(tree)
5-element Vector{Int64}:
 1
 2
 3
 4
 5
julia> treeplot(tree,a_path,porder)
347×179 Array{RGBA{N0f8},2} with eltype ColorTypes.RGBA{FixedPointNumbers.N0f8}:
 RGBA{N0f8}(1.0,1.0,1.0,1.0) ...
 ... 
```
![elimination tree of M with postorder](/img/porderpippo.png)

Then compute the row counts:
```julia
julia> rowcolcount!(M,tree,porder)
5-element Vector{Int64}:
 3
 4
 3
 2
 1
```
And even detect supernodes at the same time:
```julia
julia> rowcolcount!(M,tree,porder,true)
5-element Vector{Int64}:
 3
 4
 3
 0
 0
julia> tree
5-element Vector{Int64}:
  3
  3
  0
 -3
 -3
```
Now that the supernodal structure is reflected by the tree, plot it with [`supertreeplot`](/src/treeplots.jl) function:
```julia
julia> super_treeplot(tree,a_path,porder)
157×179 Array{RGBA{N0f8},2} with eltype ColorTypes.RGBA{FixedPointNumbers.N0f8}:
 RGBA{N0f8}(1.0,1.0,1.0,1.0) ...
 ...
```
![elimination tree of M with supernodes](/img/superpippo.png)

Get the size of supernodes :
```julia
julia> supernodecount(tree)
5-element Vector{Int64}:
 1
 1
 3
 0
 0
```

Do this all at once with the [`symbolicfact`](/src/SymbolicFact.jl) function:
```julia
julia> symbolicfact(M)
5-element Vector{Snode}:
 Snode(3, 3, -1, Set(Any[4, 3, 1]), Set{Any}(), 1)
 Snode(3, 4, -1, Set(Any[5, 4, 2, 3]), Set{Any}(), 1)
 Snode(0, 3, -1, Set(Any[5, 4, 3]), Set{Any}(), 3)
 Snode(-3, 0, -1, Set{Any}(), Set{Any}(), 0)
 Snode(-3, 0, -1, Set{Any}(), Set{Any}(), 0)
```
