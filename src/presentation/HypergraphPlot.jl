using GLMakie
using GraphMakie
using GeometryBasics
using Combinatorics


@recipe(HypergraphPlot) do scene
    Attributes(
        S_color = :yellowgreen,
        I_color = :firebrick2,
        triangle_color = (:orange, 0.3),
        hyperedge_colormap = :seaborn_bright
    )
end

function Makie.plot!(hgplot::HypergraphPlot)
    network = hgplot[1]

    # a vector of states expressed as integers 
    states = Observable(Int64[])
    # the two-section graph of the hypergraph, used as a skeleton to draw the full hypergraph
    simple_graph = Observable(Graphs.SimpleGraphs.SimpleGraph())
    # labels of the nodes
    labels = Observable(String[])
    # positions of the nodes in the image
    node_pos = Observable(Point{2, Float32}[])
    # indices of the nodes which belong to the same triangle
    faces = Observable(Matrix{Int64}(undef, 0, 3))
    
    # Called whenever the network changes to update all observables with new values
    # TODO: why is it called twice?
    function update_plot(network::HyperNetwork)
        simple_graph[] = get_twosection_graph(network)
        
        empty!(states[])
        empty!(labels[])
        faces[] = Matrix{Int64}(undef, 0, 3)
        for node = 1:get_num_nodes(network)
            push!(states[], Int64(get_state(network, node)))
            push!(labels[], "#$node, $(get_state(network, node))")
        end
        for h = 1:get_num_hyperedges(network)
            hsize = get_hyperedge_size(network, h)
            if hsize > 2
                nodes = get_nodes(network, h)
                # The nodes will be dublicated for each possible size of 
                # the hyperedge - otherwise, it is not possible to draw 
                # hyperedges of different sizes in different colors. 
                # This is why we need to offset the index of the nodes. 
                nodes .+= (hsize - 3) * get_num_nodes(network)
                # add triangular faces between all triples of nodes
                simplex_faces = collect(combinations(nodes, 3))
                simplex_faces = reduce(vcat, transpose.(simplex_faces))
                faces[] = [faces[]; simplex_faces]
            end
        end
        # trigger the update of observables
        labels[] = labels[]
        # hack to fix weird error when the whole vector contains 
        # the same value
        push!(states[], 1 - states[][end])
        states[] = states[]
        pop!(states[])
    end

    # call update_plot whenever the network changes
    Makie.Observables.on(update_plot, network)

    # call the function for the first time
    update_plot(network[])

    # map states to colors
    node_colormap = Observable{Any}()
    map!(node_colormap, hgplot.S_color, hgplot.I_color) do I_color, S_color
        [I_color, S_color]
    end

    # draw the skeleton simple graph
    gp = graphplot!(hgplot, 
                    simple_graph; 
                    nlabels=labels,
                    nlabels_distance=10,
                    node_attr = (color=states, colormap=node_colormap, markersize=12),
                    nlabels_attr = (textsize = 12, color = :gray))
    
    # graphplot gives us the positions of the nodes
    node_pos = gp[:node_pos]

    # plot the hyperedges as triangles
    # colors = [1, 1, 1, ..., 2, 2, 2, ..., 3, 3, 3, ...]
    colors = @lift repeat(1:get_max_hyperedge_size($network), inner = get_num_nodes($network))
    # stacked_node_pos = [1, 2, 3, ..., 1, 2, 3, ..., 1, 2, 3, ...]
    stacked_node_pos = @lift repeat($node_pos, outer = get_max_hyperedge_size($network))
    mesh!(hgplot,
          stacked_node_pos,
          faces,
          color = colors,
          colormap = (hgplot.hyperedge_colormap, 0.3),
          shading = false,
          transparency = true)

    return hgplot
end