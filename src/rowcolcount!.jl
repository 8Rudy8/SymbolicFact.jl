function setfind(setpath,p_leaf)
	q = p_leaf
	while (setpath[q] != q)
		q = setpath[q]
	end
	#path compression
	c = p_leaf
	while (c != q)
		tmp =setpath[c]
		setpath[c] = q
		c = tmp
	end
	return q
end

function setunion!(setpath,j,pj)
	if (pj != 0)
		setpath[j] = pj
	end
	return
end

function unflip(p)
	if (p < 0)
		unflip = -p
	else
		unflip = p
	end
	return unflip
end

function flip(p)
	if (p > 0)
		flip = -p
	else
		flip = p
	end
	return flip
end
"""
	rowcolcount!(A,parent,porder[,supernodes = false])

Return the **row count of the R factor** of `A`. If A is non-symmetric, return a `Tuple` containing the **row count** of **R** factor and the **column count** of **Q** factor.\n

## Arguments
- `A::SparseMatrixCSC{Int,Int}`: The sparse matrix on wich the count are computed.
- `parent::Vector{Int}`: The parent list of the tree.
- `porder::Vector{Int}`: The post-order of the tree.
- `supernodes::Bool`: If `true`, detect supernodes of the tree and modify the ```porder``` and ```parent``` structure to reflect the supernodes. This structure is detailed under. 

## Supernodal output structure
`porder` output will be the same as the input one, but every principal variables comes before the corresponding subordinates.
`parent` output will be on the following format :
- `parent[i] == j > 0` means that i is the principal variable of a node and j is the princpal variable of its father node. The principals variables considered here are the firsts variables of a supernode (leaves).
- `parent[i] == 0` means that i is the principal variable of a root node.
- `parent[i] == j < 0` means that i is a subordinate variable inside a node whose principal variable is `-j`.
## Example :

	         +---+
	         |7  |
	         | 6 |           
	         |  5|           parent=[7, -1, 7, -3, 0, -5, -5]
	         +---+
	        /     \\
	       /       \\
	   +--+         +--+
	   |2 |         |4 |
	   | 1|         | 3|
	   +--+         +--+

### note
The supernodes detected in this function are only linear supernodes. To detect more complexe supernodes,
an amalgamation tree is needed, wich is not implemented in this package
### See also
- [`etree`](@ref)
- [`postorder`](@ref)
"""
function rowcolcount!(A,parent,porder,supernodes::Bool=false)
	n = length(porder)
	#compute inverse post order
	iporder = zeros(Int,n)
	for j = 1:n
		iporder[porder[j]] = j
	end
	rc_R = zeros(Int,n)
	#will contain the first descendant of nodes according to their postorder
	fst_desc = fill(-1,n)
	for i = 1:n
		curr = porder[i]
		fd = curr
		if (fst_desc[curr] == -1)
			#this is a leaf node. initialize its row_count to 1
			rc_R[curr] = 1
		end
		#now we start going up the tree seeting up fd as first
		# descendant of all its ancestors
		while true
			if (fst_desc[curr] > 0)
				break
			end
			fst_desc[curr] = fd
			if (parent[curr] == 0)
				break
			end
			curr = parent[curr]
		end
	end
	if !(issymmetric(A))
		cc_Q = zeros(Int,n)
		#compute level(j), the distance from j to the root of j's subtree, for j = 1:n
		level = zeros(Int,n) 
		roots = zeros(Int,n) #this Vector will be used later to create r
		for j = 1:n
			head = parent[j]
			roots[j] = j
			while head !=0
				roots[j] = head
				level[j] = level[j]+1
				head = parent[head]
			end
		end
		#Assumme fi[i] = j means that the first nonzero row i has col index j.
		#Computing the "higher adjacency sets" (hadj).
		#The hadj of col j is defined as the union of the col indices
		#of all the rows i for wich fi(i)=j.
		#Computing in the same loops first[j], the sets
		#of rows indices where fi[i] = j for j in 1:n
		#and fi[i], the col index of the first nonzero in row i for i in 1:m
		#indlist will be used has pointers for first's subsets
		indlist = zeros(Int, n+1)
		indlist[1] = 1
		last = zeros(Int,n)
		rcnt = zeros(Int,n)
		fi = zeros(Int,A.m)
		hptr = zeros(Int,n+1)
		#we make first pass to count the number of elements in each hadj
		#and first subset, in the same loop we also compute fi
		a = 0
		for j= 1:n
			#its easier to not use porder here,
			#hadj and first does not depends of it
			add = false
			indlist[j+1] = indlist[j]
			for ii = A.colptr[j]:A.colptr[j+1]-1
				i = A.rowval[ii] #the row index
				f= fi[i]
				a = a+1
				if (f == 0) 
					add = true
					indlist[j+1] = indlist[j+1] + 1
					fi[i] = j
				elseif (j > last[f])
					rcnt[f] = rcnt[f]+1
					last[f] = j
				end	
			end
			if !add
				indlist[j+1] = indlist[j+1]+1
			end
		end
		indlist[end] = indlist[end-1]
		#transform the counts into pointers
		hptr[1] = 1
		for k = 1:n
			hptr[k+1] = hptr[k] + rcnt[k]
		end
		first = zeros(Int,indlist[end])
		hadj = zeros(Int,hptr[n+1])
		last = zeros(Int,n)
		rcnt = zeros(Int,n)
		#second loop, to fill hadj and first subsets
		for k = 1:n
			start = 1
			j = porder[k]
			#filling hadj
			for ii = A.colptr[j]:A.colptr[j+1]-1
				i = A.rowval[ii]
				f = fi[i]
				if ((k > iporder[f]) && (k > last[f]))
					ptr = hptr[f] + rcnt[f]
					hadj[ptr] = j
					rcnt[f] = rcnt[f]+1
					last[f] = k
				end
			end
			#fillling first
			for ind = indlist[j]:indlist[j+1]-1
				for i in start:A.m	
					if (fi[i]==j)
						first[ind] = i
						start = i+1
						break
					end
				end
			end
		end	
		#at this point we have everything we need to compute the rowcount of R
		#compute r[i] the smaller of i or the root of the subtree containing f_i
		r = zeros(Int,A.m)
		for i = 1:A.m
			r[i] = min(i,roots[fi[i]])
		end
		#at this point, we have everything to compute the colcount of Q
	else 
		#at this point we have everything we need to compute the rowcount
		hptr = A.colptr
		hadj = A.rowval
	end
	#setpath is used for doing path compression
	setpath = zeros(Int,n)
	for j = 1:n
		setpath[j] = j
	end
	prev_f = zeros(Int,n)
	prev_nbr = zeros(Int,n)
	#during the symbolic facto we flag as negative all the first
	#variables of a supernode. At the end, we make a simple loop
	#to modify the tree that it reflects the supernodal structure
	for jidx = 1:n
		j = unflip(porder[jidx])
		if (parent[j] != 0)
			rc_R[parent[j]] = rc_R[parent[j]]-1
		elseif(j != n)
			#if porder[j] is root, porder[j+1] is the beginning of a supernode
			porder[j+1] = flip(porder[j+1])
		end
		
		for iidx = hptr[j]:hptr[j+1]-1
			i =hadj[iidx]
			if (iporder[i] <= jidx) 
				continue
			end
			if (prev_nbr[i] == 0)
				ref = 0
			else
				ref = iporder[prev_nbr[i]]
			end
			if (iporder[fst_desc[j]] > ref)
				rc_R[j] = rc_R[j]+1
				p_leaf = prev_f[i]
				porder[jidx] = flip(porder[jidx])
				if (p_leaf != 0)
					q = setfind(setpath,p_leaf)
					#q is a leaf or a row subtree, thus the first variable in 
					#a supernode. flag it
					porder[iporder[q]] = flip(porder[iporder[q]])
					rc_R[q] = rc_R[q] - 1
					#u = porder[i] ?
					#ou u = i
				end
				prev_f[i] = j
			end
			prev_nbr[i] = j
		end
		setunion!(setpath,j,parent[j])
		if !issymmetric(A)
			for ind = indlist[j] : indlist[j+1]-1
				i = first[ind]
				if (i ==0)
					continue
				end
				cc_Q[j] = cc_Q[j] +1
				if (parent[r[i]] != 0) 
					#if r[i] is not the root of the subtree
					cc_Q[parent[r[i]]] = cc_Q[parent[r[i]]] - 1
				end
			end
		end
	end
	for jj= 1:n-1
		j = unflip(porder[jj])
		if (parent[j] != 0)
			rc_R[parent[j]] = rc_R[parent[j]] + rc_R[j]
			if !issymmetric(A)
				cc_Q[parent[j]] = cc_Q[parent[j]] + cc_Q[j]
			end
		end
	end
	
	if supernodes
		#supernode detection
		#now we have to update the tree in order to reflect the supernodal structure.
		#we simply go through the postorder looking for the variables we flagged above
		i = 1
		while (i < n)
			porder[i] = unflip(porder[i])
			j = i+1
			while ((j <= n) && (porder[j] > 0))
				j = j+1
			end
			#at this point, j is the first variable of the next supernode and
			#j-1 is the principal variable of the supernode starting at i.
			#all the variables i:j-2 are subordinate to j-1
			parent[porder[i]] = parent[porder[j-1]]
			for k = (i+1):(j-1)
				parent[porder[k]] = -porder[i]
				rc_R[porder[k]] = 0	
			end
			#we want the supervariable to be in first position in the porder
			i = j
		end
		porder[n] = unflip(porder[n])
	else
		#unflag the porder
		for i = 1:n
			porder[i] = unflip(porder[i])
		end
	end
	if !issymmetric(A)
		return rc_R,cc_Q
	else
		return rc_R
	end
end
