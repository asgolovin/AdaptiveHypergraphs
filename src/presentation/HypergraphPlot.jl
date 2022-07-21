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
    colors = Observable(Int64[])
    # dummy example with scatter 
    xs = Observable(Int64[])
    ys = Observable(Int64[])

    function update_plot(network::HyperNetwork)
        println("Plot has changed!")
        # update the observables
        xs[] = 1:get_num_nodes(network)
        ys[] = [get_node_degree(network, i) for i in xs[]]
        
        empty!(colors[])
        for node = 1:get_num_nodes(network)
            push!(colors[], Int(get_state(network, node)))
        end
        colors[] = colors[]
    end

    # call update_plot whenever the network changes
    Makie.Observables.on(update_plot, network)

    # call the function for the first time
    update_plot(network[])

    colormap = Observable{Any}()
    map!(colormap, hgplot.S_color, hgplot.I_color) do I_color, S_color
        [I_color, S_color]
    end

    # dummy example with scatter - replace later!
    scatter!(hgplot, xs, ys; color=colors, colormap=colormap)

    return hgplot
end

n = 10
network = Observable(HyperNetwork(n, 0.3))
build_RSC_hg!(network[], (20, 4))

f = Figure()
display(f)

hypergraphplot(f[1, 1], network)