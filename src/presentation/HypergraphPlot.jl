using GLMakie

include("../data/Network.jl")

@recipe(HypergraphPlot) do scene
    Attributes(
        S_color = :green,
        I_color = :red
    )
end

function Makie.plot!(hgplot::HypergraphPlot)
    network = hgplot[1]

    # define plot elements that can change with time as observables
    # a vector of states expressed as integers 
    states = Observable(Int64[])
    # the two-section graph of the hypergraph, used as a skeleton to draw the full hypergraph
    simple_graph = Observable(Graphs.SimpleGraphs.SimpleGraph())
    
    # Called whenever the network changes to update all observables with new values
    function update_plot(network::HyperNetwork)
        println("Plot has changed!")
        simple_graph[] = get_twosection_graph(network)
        
        empty!(states[])
        for node = 1:get_num_nodes(network)
            push!(states[], Int(get_state(network, node)))
        end
        states[] = states[]
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
    graphplot!(hgplot, simple_graph)

    return hgplot
end

n = 10
network = Observable(HyperNetwork(n, 0.3))
build_RSC_hg!(network[], (20, 4))

f = Figure()
display(f)

hypergraphplot(f[1, 1], network)