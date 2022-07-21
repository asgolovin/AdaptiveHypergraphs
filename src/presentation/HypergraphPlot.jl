using GLMakie

include("../data/Network.jl")

@recipe(HypergraphPlot) do scene
    Attributes(
        S_color = :green,
        I_color = :red,
        hyperedge_colormap = :GnBu
    )
end

function Makie.plot!(hgplot::HypergraphPlot)
    network = hgplot[1]

    # define plot elements that can change with time as observables

    # a vector of states expressed as integers 
    states = Observable(Int64[])
    # the two-section graph of the hypergraph, used as a skeleton to draw the full hypergraph
    simple_graph = Observable(Graphs.SimpleGraphs.SimpleGraph())
    # labels of the nodes
    labels = Observable(String[])
    # positions of the nodes in the image
    node_pos = Observable(Point{2, Float32}[])
    # indices of the nodes which belong to the same triangle
    faces = Observable(Int64[])
    
    # Called whenever the network changes to update all observables with new values
    # TODO: why is it called twice?
    function update_plot(network::HyperNetwork)
        simple_graph[] = get_twosection_graph(network)
        
        empty!(states[])
        empty!(labels[])
        faces[] = []
        for node = 1:get_num_nodes(network)
            push!(states[], Int(get_state(network, node)))
            push!(labels[], "#$node, $(get_state(network, node))")
        end
        for h = 1:get_num_hyperedges(network)
            if get_hyperedge_size(network, h) > 2
                nodes = get_nodes(network, h)
                faces[] = [faces[]; nodes]
            end
        end
        states[] = states[]
        labels[] = labels[]
    end

    # call update_plot whenever the network changes
    Makie.Observables.on(update_plot, network)

    # call the function for the first time
    update_plot(network[])

    # map states to colors
    colormap = Observable{Any}()
    map!(colormap, hgplot.S_color, hgplot.I_color) do I_color, S_color
        [I_color, S_color]
    end

    # draw the skeleton simple graph
    gp = graphplot!(hgplot, 
                    simple_graph; 
                    nlabels=labels,
                    nlabels_distance=10,
                    node_attr = (color=states, colormap=colormap, markersize=15))
    
    # graphplot gives us the positions of the nodes
    node_pos = gp[:node_pos]

    # plot the hyperedges as triangles
    mesh!(hgplot, node_pos, faces, color = (:orange, 0.3), shading = false)

    return hgplot
end

n = 10
network = Observable(HyperNetwork(n, 0.3))
build_RSC_hg!(network[], (20, 5))

f = Figure()
display(f)

hypergraphplot(f[1, 1], network)