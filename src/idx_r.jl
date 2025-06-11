"""
	isleaf(node::Int,parent::Vector{Int})

**Return** a `Bool`, `true` if `node` is a leaf in the tree `parent`, else `false`.
"""
function isleaf(node,parent)
	isleaf= true
	for i = 1:length(parent)
		if (parent[i] == node)
			#node has a son, he is not a leaf
			isleaf = false
		end
	end
	if (parent[node] <= 0)
		#node not the first in a supernode, he is not a leaf neither
		isleaf = false
	end
	return isleaf
end
"""
	idx_r(A,p,porder)

Return the row indices of each non-zeros element per column, fill-in coef included.

## Arguments

- `A::SparseMatrixCSC{Float64,Int}`: The matrix computed, have to be structurally symmetric.
- `p::Vector{Int}`: The elimination tree of the matrix, can reflect or not supernodes.
- `porder::Vector{Int}`: The post-order of the tree.

### Note
For reasons of optimization, this function **does not work for symmetric matrices** yet.

### See also
- [`etree`](@ref)
- [`porder`](@ref)
"""
function idx_r(A,p,porder)
	@assert issymmetric(A)
	n = length(p)
	prev= []
	I=[Set() for i=1:n]
	update=[Set() for i=1:n]
	for k =1:n
		j = porder[k]
		if isleaf(j,p)
			#the node is a leaf, they will not be filled-in
			#we can directly compute their row indices
			I[j] = Set(i for i in A[:,j].nzind if (i >= j))
		end
		if (p[j] < 0)
			#this node is a part of a supernode
			continue
		else
			#comupting row indices : the union of the row indices of j
			#and the row indices of the sons of j.
			I[j] = union(Set(i for i in A[:,j].nzind if (i >= j)),update[j])
		end
		if (p[j] != 0)
			#keeping the rows to add to the parent node in the future
			update[p[j]] = intersect(Set(i for i=p[j]:n),union(update[p[j]],I[j]))
		end
	end
	return I
end
