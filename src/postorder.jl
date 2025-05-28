"""
		postorder(parent)

**Returns** an array that gives the  post-order of the nodes of the `parent` tree, named ```porder```, 
defined by the following property : 

- ```porder[k] == i``` means that i is the k-th node in the post-order.
### See also
- [`etree`](@ref)
"""
function postorder(parent)

	n = length(parent)
	son = zeros(Int,n)
	brother = zeros(Int,n)
	stack = zeros(Int,n)
	porder = zeros(Int,n)
	#construction of a brother and son Vector
		for i = n:-1:1
			father = parent[i]
			if (father != 0)
				br = son[father]
				brother[i] = br
				son[father] = i
			end
		end
	#dfs traversing
	head = 0
	hp = 0
	pp = 1
	for t = 1:n
		if (parent[t] != 0)
			continue
		end 
		#in this case, t is the root of a tree, will be executed for each tree of the forest
		hp = hp+1
		stack[hp] = t	
		head = t
		while true
			if (son[head]==0)
				#we reached the bottom
				porder[pp] = head
				pp = pp+1
				hp = hp-1
				if (parent[head] != 0)
					son[parent[head]] = brother[head]
				end
				if (hp == 0)
					#tree is exhausted
					break
				end
				head = stack[hp]
			else
				#go down one more level
				hp = hp + 1
				stack[hp] = son[head]
				head = son[head]
			end
		end
	end
	return porder
end
