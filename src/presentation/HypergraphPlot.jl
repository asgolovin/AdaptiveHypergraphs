using GLMakie
using GeometryBasics
using Combinatorics

include("../data/Network.jl")

@recipe(HypergraphPlot) do scene
    Attributes(
        S_color = :yellowgreen,
        I_color = :firebrick2,
        triangle_color = (:orange, 0.3),
        hyperedge_colormap = :GnBu
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
    # for each face stores the size of the hyperedge to which the face belongs to
    hyperedge_size = Observable(Int64[])
    
    # Called whenever the network changes to update all observables with new values
    # TODO: why is it called twice?
    function update_plot(network::HyperNetwork)
        simple_graph[] = get_twosection_graph(network)
        
        empty!(states[])
        empty!(labels[])
        faces[] = Matrix{Int64}(undef, 0, 3)
        empty!(hyperedge_size[])
        for node = 1:get_num_nodes(network)
            push!(states[], Int(get_state(network, node)))
            push!(labels[], "#$node")
        end
        for h = 1:get_num_hyperedges(network)
            hsize = get_hyperedge_size(network, h)
            if hsize > 2
                nodes = get_nodes(network, h)
                # add triangular faces between all triples of nodes
                simplex_faces = collect(combinations(nodes, 3))
                simplex_faces = reduce(vcat, transpose.(simplex_faces))
                faces[] = [faces[]; simplex_faces]
                push!(hyperedge_size[], hsize)
            end
        end
        # trigger the update of observables
        states[] = states[]
        labels[] = labels[]
        hyperedge_size[] = hyperedge_size[]
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
                    node_attr = (color=states, colormap=node_colormap, markersize=15),
                    nlabels_attr = (textsize = 12, color = :gray))
    
    # graphplot gives us the positions of the nodes
    node_pos = gp[:node_pos]

    # plot the hyperedges as triangles
    mesh!(hgplot, 
          node_pos,
          faces,
          color = (:orange, 0.3),
          #colormap = hgplot.hyperedge_colormap,
          shading = false)

    return hgplot
end

n = 10
network = Observable(HyperNetwork(n, 0.3))
build_RSC_hg!(network[], (10, 2, 1))

f = Figure()
display(f)

hypergraphplot(f[1, 1], network)