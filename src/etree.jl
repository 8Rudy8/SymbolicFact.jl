    """
    add_node!(ancestor,parent,j,i)
        
Adds `j` as the root of the path that contains `i` in the `parent` tree, if it is not already in it.
    """
function add_node!(ancestor, parent, j, i)
    if(i==0)
        # i is the root of the tree, do nothing
        return
    end
    k = i
    while true
    r = ancestor[k]
    if (r==j)
        # node was already added, do nothing
        return
    end
    ancestor[k] = j
    if(r == 0)
        # we reached the top of the subtree
        # add j on top
        parent[k] = j
        return
    end
    k = r
    end
end
"""
etree(A)

This function returns the elimination tree of the matrix A, `parent`, defined by the following properties :

- `parent(i) == j`  means that `j` is the parent of  `i`.

- `parent(i) == 0`, means that `i` has no parents.

## Arguments
- `A::SparsematrixCSC{Int,Int}`: the matrix used to create the elimination Tree. If `A` is non-symmetric, the tree computed will be the one of produce of  transposed `A`  and `A`.
"""
function etree(A)
    #ancestor is used for doing path compression
    ancestor = zeros(Int, A.n)
    parent   = zeros(Int, A.n)
    if(!issymmetric(A))
        prev_col = zeros(Int, A.m)
    end
    for j = 1:A.n
        for iidx = A.colptr[j]:A.colptr[j+1]-1
            i = A.rowval[iidx]
            if(issymmetric(A))
                if (i >= j)
                    continue
                end
                k = i
            else
                k = prev_col[i]
            end
            add_node!(ancestor, parent, j, k)
            if(!issymmetric(A)) 
                prev_col[i] = j
            end
        end
    end
    return parent
end