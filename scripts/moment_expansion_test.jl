using AdaptiveHypergraphs
using GLMakie
using ColorSchemes

include("../input/large_network.jl")
nparams = params.network_params
mparams = params.model_params
bparams = params.batch_params
vparams = params.visualization_params

# SOLVING

n = nparams.num_nodes
max_size = length(nparams.num_hyperedges) + 1
network = HyperNetwork(n, nparams.infected_prob, max_size)
build_RSC_hg!(network, nparams.num_hyperedges)
motif_count = get_motif_count(network)

t, sol = moment_expansion(motif_count, params, (0., 1000.))

# PLOTTING

fig = Figure()
ax = Axis(fig[1, 1])

node_colormap = params.visualization_params.node_colormap

linecolors = Dict()
linestyles = Dict()
colorscheme = colorschemes[node_colormap]
for label in keys(sol)
    numA = label.left_total[AdaptiveHypergraphs.A] + label.right[AdaptiveHypergraphs.A]
    numB = label.left_total[AdaptiveHypergraphs.B] + label.right[AdaptiveHypergraphs.B]
    ratio = numB / (numA + numB)
    linecolor = get(colorscheme, ratio)
    linecolors[label] = linecolor
    linestyles[label] = [nothing, :dash, :dot, :dashdot, :dashdotdot][numA + numB]
end

for label in keys(sol)
    lines!(t, sol[label], color=linecolors[label], label = "$label", linestyle = linestyles[label])
end

Legend(fig[1, 2], ax)

fig