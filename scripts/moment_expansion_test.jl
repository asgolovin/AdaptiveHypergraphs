using AdaptiveHypergraphs
using GLMakie
using ColorSchemes
using DifferentialEquations

include("../input/large_network.jl")

Label = AdaptiveHypergraphs.Label

# rhs only for the equations for a simple graph, not a hypergraph
function rhs_to_random!(dx, x, p, t)
    A, B, AA, BB, AB = x

    dx[1] = 0.0
    dx[2] = 0.0

    # closures
    ABA = 0.5 * AB^2 / B
    ABB = 2.0 * BB * AB / B
    AAB = 2.0 * AA * AB / A
    BAB = 0.5 * AB^2 / A

    dx[3] = 0.5 * AB * (1 - p + p * A / n) + 0.5 * (1 - p) * (2.0 * ABA - AAB)
    dx[4] = 0.5 * AB * (1 - p + p * B / n) + 0.5 * (1 - p) * (2.0 * BAB - ABB)
    dx[5] = -dx[3] - dx[4]

    return dx
end

# rhs only for the equations for a simple graph, not a hypergraph
function rhs_to_same!(dx, x, p, t)
    A, B, AA, BB, AB = x

    dx[1] = 0.0
    dx[2] = 0.0

    # closures
    ABA = 0.5 * AB^2 / B
    ABB = 2.0 * BB * AB / B
    AAB = 2.0 * AA * AB / A
    BAB = 0.5 * AB^2 / A

    dx[3] = 0.5 * AB + 0.5 * (1 - p) * (2.0 * ABA - AAB)
    dx[4] = 0.5 * AB + 0.5 * (1 - p) * (2.0 * BAB - ABB)
    dx[5] = -dx[3] - dx[4]

    return dx
end

nparams = params.network_params
mparams = params.model_params
bparams = params.batch_params
vparams = params.visualization_params

# SOLVING

n = nparams.num_nodes
max_size = length(nparams.num_hyperedges) + 1
network = HyperNetwork(n, nparams.infected_prob, max_size)
build_RSC_hg!(network, nparams.num_hyperedges)
tspan = (0.0, 1000.0)

t, sol = moment_expansion(params, tspan, moment_closure)

if max_size == 2
    x0 = [motif_count[Label("[A]")],
          motif_count[Label("[B]")],
          motif_count[Label("[A2]")],
          motif_count[Label("[B2]")],
          motif_count[Label("[AB]")]]
    p = mparams.adaptivity_prob
    adaptivity_rule = mparams.adaptivity_rule

    if typeof(adaptivity_rule) <: RewireToRandom
        problem = ODEProblem(rhs_to_random!, x0, tspan, p)
    elseif typeof(adaptivity_rule) <: RewireToSame
        problem = ODEProblem(rhs_to_same!, x0, tspan, p)
    end
    sol_simple = solve(problem)
    t_simple = sol_simple.t
    u_simple = reduce(hcat, sol_simple.u)
end

# PLOTTING

fig2 = Figure()
ax = Axis(fig2[1, 1])

node_colormap = params.visualization_params.node_colormap

linecolors = []
linestyles = []
colorscheme = colorschemes[node_colormap]
for (i, label) in enumerate(keys(sol))
    AB_ratio = label.left_total[AdaptiveHypergraphs.A] / (size(label)[1])
    push!(linecolors, get(colorscheme, AB_ratio, (0.0, 1.0)))
    push!(linestyles, [nothing, :dash, :dot][size(label)[1]])
end

for (i, label) in enumerate(keys(sol))
    lines!(t, sol[label]; label="$label", color=linecolors[i], linestyle=linestyles[i],
           linewidth=3.0)
end

if max_size == 2
    for (i, label) in enumerate(["A", "B", "AA", "BB", "AB"])
        lines!(t_simple, u_simple[i, :]; label=label, color=linecolors[i])
    end
end

ylims!(-100.0, nothing)

Legend(fig2[1, 2], ax)

println("AB: $(sol[AdaptiveHypergraphs.Label("[AB]")][end])")
println("AA: $(sol[AdaptiveHypergraphs.Label("[A2]")][end])")
println("BB: $(sol[AdaptiveHypergraphs.Label("[B2]")][end])")

fig2