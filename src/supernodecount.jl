"""
	supernodecount(p)

Count the number of node in each supernode of `p`, **return** the result with a list of counts per nodes :

- `supernodecount(p)[i] == j` means that the supernode `i` is made of `j` nodes.
- `supernodecount(p)[i] == 0` means that the node `i` is a subnode.
### See also
- [`etree`](@ref)
- [`rowcolcount!`](@ref)
"""
function supernodecount(p)
	n = length(p)
	np = zeros(Int,n)
	for i=1:n
		if p[i]<0
			#i is the subnode of the supernode -p[i]
			np[-p[i]] = np[-p[i]]+1
		end
	end
	for i=1:n
		if (p[i] >= 0)
			np[i] = np[i] + 1
		end
	end
	return np 
end
