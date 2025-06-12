module SymbolicFact
	using MatrixDepot
	using LinearAlgebra
	using SparseArrays
	using Plots
	using GraphvizDotLang: graph, digraph, node, edge, attr, save, subgraph
	using Metis
	using Images, FileIO
	using D3Trees
    
    export 
        #functions
        etree,
        fill_reducing_perm,
        postorder,
        rowcolcount!,
        idx_r,
		symbolicfact,
		supernodecount,
		#plotting functions
		treeplot,
		super_treeplot,
		D3_treeplot,
        #types
        Snode

    include("etree.jl")
    include("fill_reducing_perm.jl")
    include("postorder.jl")
    include("rowcolcount!.jl")
    include("idx_r.jl")
	include("supernodecount.jl")
	include("treeplots.jl")
"""
	Snode(father::Int,rc::Int,cc::Int,idxr::Set,idxc::Set,np::Int,porder::Int)
`DataType` made for containing the informations of symbolic factorization. Represent a node or a supernode of an elimination tree.

## Attributes
- `father::Int`: The father of the considered node in the elimmination tree, can reflect or not the supernodal structure accordeing to the rules described in [`rowcolcount!`](@link) doc.
- `rc::Int`: The **row count** of the considered node, equal `0` if the node is a subnode.
- `cc:Int`: The **column count** of the considered node, equal `-1` if there is no column count.
- `idxr::Set`: The **row indices** of the left column of the front matrix. Equal `Set([])` if the node considered is a subnode.
- `idxc::Set`: Not implemented
- `np::Int`: the supernode size, equal 0 if the node is a subnode
- `porder::Int`: the post-order of the node

### See also
- [`etree`](@link)
- [`rowcolcount!`](@link)
"""
struct Snode
	father ::Int 
	rc::Int
	cc::Int
	idxr::Set
	idxc::Set
	np::Int
	porder::Int
	#defining constructors
	Snode(p,x,y,z,t,e,a) = new(p,x,y,z,t,e,a) #with all params
	Snode() = new(0,0,0,Set(),Set(),1,0)#empty node
end
"""
	symbolicfact(A[,supernodes = true])

Make the symbolic factorization of `A`.

## Arguments
- `A::SparseMatrixCSC{Float64,Int}`: the sparse matrix that will be factorized
- `supernodes::Bool`: If `true`, the structure returned will reflect supernode detected by the algo. See [`rowcolcount!`](@link) for the detail.

## Return
Return a `Vector` containing the symbolic factorization informations. Each element of it is an instance of [`Snode`](@link) representing a column.

### See also
- [`Snode`](@link)
- [`rowcolcount!`](@link)
"""
function symbolicfact(A,supernodes = true)
	elimtree = etree(A)
	porder = postorder(elimtree)
	n = length(porder)
	result = [Snode() for i=1:n]
	if issymmetric(A)
		rowcount = rowcolcount!(A,elimtree,porder,supernodes)
		np = supernodecount(elimtree)
		idxr = idx_r(A,elimtree,porder)
		for i = 1:n
			#we set the column count to -1 to make difference beetween
			#no column count and column count equal 0
			result[i] = Snode(elimtree[i],rowcount[i],-1,idxr[i],Set(),np[i],findall(x -> x== i,porder)[1])
		end
	else
		rowcount,colcount = rowcolcount!(A,elimtree,porder,supernodes)
		np = supernodecount(elimtree)
		#idxr = idx_r(B,elimtree,porder) #not implemented yet for a non-symmetric matrix
		for i = 1:n
			result[i] = Snode(elimtree[i],rowcount[i],colcount[i], Set(),Set(),np[i],findall(x -> x== i,porder)[1])
		end
	end
	return result
end
end # module SymbolicFact
