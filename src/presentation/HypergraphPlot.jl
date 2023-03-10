using GLMakie
using GraphMakie
using GeometryBasics
using Combinatorics
using ColorSchemes

export hypergraphplot, hypergraphplot!

"""
    hypergraphplot(network::HyperNetwork)
    hypergraphplot!(ax, network::HyperNetwork)

Plots a visualization of a hypergraph with nodes in particular states. 

The visualization consists of three consecutive steps:

    1. First, the two-section of the hypergraph is plotted using the `graphplot()` function 
    from GraphMakie. A two-section of a hypergarph is a simple graph in which two nodes 
    are connected if they belong to the same hyperedge. Importantly, the `graphplot()` 
    function also computes the positions of the nodes on a two-dimensional plane. Those 
    node positions are used for the two following steps. 

    2. Having obtained the node positions, hyperedges of size greater than two are plotted 
    using meshes: within each hyperedge, all triples of nodes are connected with a triangle. 
    This produces a visualization similar to what is usually used to plot simplical complexes: 
    hyperdeges with three nodes are visualized as triangles, hyperdeges with four nodes 
    as tetrahedra and so on. Hyperedges of different sizes are plotted in different colors. 

    3. Finally, hyperdeges with two nodes, i.e., simple edges, are plotted as lines. The 
    color of the lines is sampled from the hyperdege colormap. 
"""
GLMakie.@recipe(HypergraphPlot) do scene
    return Attributes(; node_colormap=:RdYlGn_6,
                      hyperedge_colormap=:thermal)
end

function GLMakie.plot!(hgplot::HypergraphPlot)
    network = hgplot[1]

    # a vector of node states expressed as integers 
    states = Observable(Int64[])
    # the two-section graph of the hypergraph, used as a skeleton to draw the full hypergraph
    simple_graph = Observable(Graphs.SimpleGraphs.SimpleGraph())
    # labels of the nodes
    labels = Observable(String[])
    # coordinates of the nodes in the image
    node_pos = Observable(Point{2,Float32}[])
    # indices of the nodes which belong to the same triangle
    faces = Observable(Matrix{Int64}(undef, 0, 3))
    # indices of the nodes which belong to the same edge
    edges = Observable(NTuple{2,Int64}[])
    # historical maximum hyperdege size
    max_hedge_size = get_max_size(network[])
    # coordinates of the nodes stacked on top of eack other max_hedge_size times.
    # We need this artificial dublication to draw hyperedges of different sizes in different colors.
    # stacked_node_pos = [1, 2, 3, ..., 1, 2, 3, ..., 1, 2, 3, ...]
    stacked_node_pos = Observable(Point{2,Float32}[])
    # colors of the triangle faces
    # colors = [1, 1, 1, ..., 2, 2, 2, ..., max_hedge_size, max_hedge_size, max_hedge_size, ...]
    colors = repeat(1:max_hedge_size; inner=get_num_nodes(network[]))

    # Called whenever the network changes to update observables with new values
    Makie.Observables.on(network) do network
        simple_graph[] = get_twosection_graph(network)

        empty!(states[])
        empty!(labels[])
        empty!(edges[])
        for node in get_nodes(network)
            push!(states[], Int64(get_state(network, node)))
            push!(labels[], "#$node, $(get_state(network, node))")
        end
        for h in get_hyperedges(network)
            hsize = get_hyperedge_size(network, h)
            nodes = get_nodes(network, h)
            if hsize == 2
                # add edges between the nodes
                push!(edges[], Tuple(nodes))
            end
        end
        # trigger the update of observables
        notify(labels)
        # hack to fix a weird error when the whole vector contains 
        # the same value
        push!(states[], 1 - states[][end])
        notify(states)
        pop!(states[])
        return nothing
    end

    # call the function for the first time
    notify(network)

    # draw the skeleton simple graph
    gp = graphplot!(hgplot,
                    simple_graph;
                    nlabels=labels,
                    nlabels_distance=10,
                    edge_color=:gray,
                    edge_width=1.5,
                    node_attr=(color=states, colormap=hgplot.node_colormap,
                               markersize=12),
                    nlabels_attr=(textsize=12, color=:gray))

    # graphplot gives us the positions of the nodes
    node_pos = gp[:node_pos]

    # using the positions of the nodes, update Observables related to the triangle faces
    Makie.Observables.onany(network, node_pos) do network, node_pos
        faces.val = Matrix{Int64}(undef, 0, 3)
        for h in get_hyperedges(network)
            hsize = get_hyperedge_size(network, h)
            nodes = get_nodes(network, h)
            if hsize > 2
                # The indices are offsetted because they refer to the 
                # position in the stacked_node_pos vector, not the node_pos 
                # vector. 
                nodes .+= (hsize - 2) * get_num_nodes(network)
                # add triangular faces between all triples of nodes
                simplex_faces = collect(combinations(nodes, 3))
                simplex_faces = reduce(vcat, transpose.(simplex_faces))
                faces.val = [faces[]; simplex_faces]
            end
        end
        stacked_node_pos[] = repeat(node_pos; outer=max_hedge_size)
        return nothing
    end

    notify(node_pos)

    mesh!(hgplot,
          stacked_node_pos,
          faces;
          color=colors,
          colormap=(hgplot.hyperedge_colormap, 0.5),
          shading=false,
          transparency=true)

    edge_node_pos = @lift [($node_pos[first], $node_pos[second])
                           for (first, second) in $edges]

    # plot the 2-hyperedges
    linesegments!(hgplot,
                  edge_node_pos;
                  color=colorschemes[hgplot.hyperedge_colormap[]][1])

    return hgplot
end