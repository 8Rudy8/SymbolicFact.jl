"""
	treeplot(p[,porder])

Plot the tree/forest `p` thanks to the `GraphvitzDotLang` pkg.\n
The resolution of the plotted graph is limited, it is not adviced to use
it for graph with more than **50** nodes.
## Arguments
- `p::Vector{Int}`: the parent `Vector` of the tree
- `porder::Vector{Int}`: the post-order of the tree, if specified, each node value will be plotted with its port-order number with the format '**value : number**'.
### See also
- [`etree`](@ref)
- [`postorder`](@ref)
"""
function treeplot(p,porder=[])
	#creating graph
	g = digraph(rankdir="BT")
	if (porder == [])
		for i=1:length(p)
			if (p[i] != 0)
				#adding edges
				g |> edge(string(i),string(p[i]))
			end
		end
	else	
		for i=1:size(p,1)
		#adding nodes
		g |> node(string(i); label = string(i)* " : "*string(findall(x -> x== i,porder)[1]))
			if (p[i] != 0)
				g |> node(string(p[i]))
				#adding edges
				g |> edge(string(i),string(p[i]))
			end	
		end
	end
	save(g,"/tmp/pippo.png",engine="dot", format="png")
	img = load("/tmp/pippo.png")
end

"""
	super_treeplot(p,porder[,showall = false])

Plot the elimination tree contained in `p`, with **`p` reflecting supernodes**.\n

## Arguments
- `p::Vector{Int}`: Represents parent list of the tree, can reflect supernodes with the structure given by the [`rowcolcount!`](@link) function. In this case, supernodes will be plotted with the format **first node** : **subnodes**.
- `porder::Vector{Int}`: The post-order of the tree.
- `showall::Bool`: If `true`, each values of the subnodes are plotted. If `false`, just the first and the last subnode of a supernode are plotted.
### See also
- [`etree`](@ref)
- [`postorder`](@ref)
- [`rowcolcount!`](@ref)
"""
function super_treeplot(p,porder,showall:: Bool = false)
	#creating graph
	g = digraph(rankdir="BT")
	i = 1
	n = length(porder)
	lab = ""
	for i = 1:n
		if (p[porder[i]] >= 0)
			#porder[i] is the first value of a supernode
			main = porder[i]
			name = string(main)
			others = []
			already_value = false
			#looking for others values  of the supernode porder[i]
			for j = 1:n
				if (p[porder[j]] == -main)
					if (!already_value)
						name = name*": "
						already_value = true
					end
					push!(others , string(porder[j]))
				end
			end
			#build the name
			if ((already_value) && (showall))
				name = name*join(others,", ") 
			elseif already_value
				if length(others)> 3
					name = name*string(others[begin])*", ... , "*string(others[end])
				else
					name = name*join(others,", ")
				end
			end
			#add the node and link it
			g |> node(string(main);label = name,shape = "record", style ="rounded")
			if (p[porder[i]] > 0)
				g |> edge(string(main), string(p[main]))
			end
		end
	end
	#plot the graph
	save(g,"/tmp/pippo.png",engine="dot", format="png")
	img = load("/tmp/pippo.png")
end

"""
	D3_treeplot(parent,nav="firefox"[;porder=[],title="Tree"])

Plot the Tree described by the `parent` list in you navigator thanks to the `D3Trees` pkg.\n 

## Arguments
- `parent::Vector{Int}`: The parent `Vector` of the tree.
- `nav::String`: the name of the navigator used to plot the tree, do not forget to precise wich navigator you are using if it's not `"firefox"`. To get more details, see [`D3Trees.jl`](@link) doc.

## Keywords arguments
- `porder::Vector{Int}`: If specified, the post-order of a node will be plotted while hovering it.
- `title::String`: The HTML title of the plotted tree.

### Note
This function is **not able to plot forests** yet.

### See also
- [`etree`](@ref)
- [`postorder`](@ref)
- [`D3Trees.jl`](@ref)
"""
function D3_treeplot(parent,nav::String="firefox";porder::Vector{Any}=[],title::String="Tree")
	##TO UPDATE FOR PLOTTING FORESTS
	n = length(parent)
	roots = [] #will contain the roots of each tree
	#compute the children Vector of the tree(s)
	children = [[] for i=1:n]
	for i = 1:n
		if (parent[i] == 0)
			push!(roots,i)
		else
			push!(children[parent[i]],i)
		end
	end

	#we can't plot the tree from a leaf to the root (the root have to be in first position)
	#we are establishing a new order for plotting the tree from the root to the leaf
	num = zeros(Int,n)
	#roots have to be in first position
	for i =1:length(roots)
		num[i] = roots[i]
	end
	heads = roots
	i = length(roots)
	while i < n
		new_heads = []
		for j in heads
			push!(new_heads,children[j])
		end
		heads = collect(Iterators.flatten(new_heads))
		for k in 1:length(heads)
			i = i+1
			num[i] = convert(Int64,heads[k])
		end
	end
	#computing a new children Vector of the tree(s) where roots are in first position
	new_children = [[] for i = 1:n]
	for i = 1:n
		k = num[i]
		for elt in children[k]
			push!(new_children[i],findall(x->x==elt,num)[1])
		end
	end

	#tree params
	click_size = 2
	if (n > 30)
		click_size = 4
	end

	if (porder == [])
		tooltip = []
	else
		tooltip = [string(findall(x->x==i,porder)[1]) for i in num]
	end
	t = D3Tree(new_children,
			   init_expand=n,
			   title = title,
			   text =[string(i) for i in num],
			   tooltip = tooltip,
			   init_duration = 500,
			   svg_height = 400,
			   on_click_display_depth= click_size)
	#plotting the tree
	if (nav == "chrome")
		inchrome(t)
	else
		inbrowser(t,nav)
	end
end
