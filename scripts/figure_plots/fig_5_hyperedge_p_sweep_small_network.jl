using StatsBase

using AdaptiveHypergraphs
const AH = AdaptiveHypergraphs

include("figure_plotting_tools.jl")

# ======================================================
# -------------------- INPUT ---------------------------

draw_error_bars = true

input_folder = joinpath(projectdir(),
                        "results/run_2023-06-16_09-30-39_p_sweep_small_network")

panels = [:maj_rts, :maj_rtr, :prop_rts, :prop_rtr]
titles = Dict(:maj_rts => "Majority voting,\nrewire-to-same",
              :maj_rtr => "Majority voting,\nrewire-to-random",
              :prop_rts => "Proportional voting\nrewire-to-same",
              :prop_rtr => "Proportional voting\nrewire-to-random")

num_nodes = 10000
max_size = 4
linecolors = hyperedge_linecolors(max_size)

prompt = false

# number of active hyperedges below which the simulation counts as converged
converged_threshold = 10

filename = "./figures/fig_5_hyperedges_p_sweep_small_network.pdf"

xlabel = L"$p$"
ylabel = L"\mathbf{M}(t)"
xticks = 0:0.2:1
yticks = 0:250:1200
xlims = (-0.02, 1.02)
ylims = (0, 1250)

max_Δp = 0.5
max_Δh = 80

# ------------------------------------------------------
# ======================================================

# create figure
fig = create_figure(:large, 3 / 2)
ax_gridpos = fig[1, 1] = GridLayout()

axes_matrix = add_four_rules_axes!(ax_gridpos, titles,
                                   xlabel, ylabel,
                                   xticks, yticks,
                                   xlims, ylims)

p_sweep = collect(0:0.05:1)
hyperedge_dist = Dict(panel => Dict(s => Float64[] for s in 2:max_size)
                      for panel in panels)
hyperedge_dist_error = Dict(panel => Dict(s => Float64[] for s in 2:max_size)
                            for panel in panels)

data_folder = DataFolder(input_folder, :simulation)
for (batch_folder, batch_num) in data_folder
    params = load_params(joinpath(batch_folder.folder, "input_params.json"))

    prop_rule = params.model_params.propagation_rule
    prop_symb = typeof(prop_rule) <: ProportionalVoting ? :prop : :maj
    adapt_rule = params.model_params.adaptivity_rule
    adapt_symb = typeof(adapt_rule) <: RewireToSame ? :rts : :rtr

    panel = Symbol("$(prop_symb)_$(adapt_symb)")

    p = params.model_params.adaptivity_prob

    final_hyperedges = Dict(s => Float64[] for s in 2:max_size)
    for (run_folder, run_num) in batch_folder
        for s in 2:max_size
            _, hyperedge_count = load_data("hyperedge_count_$s.csv", run_folder)

            # take the last 10% of the solution and average it
            num_points = length(hyperedge_count)
            start = Int64(round(num_points * 0.9))
            push!(final_hyperedges[s], mean(hyperedge_count[start:end]))
        end
    end

    for s in 2:max_size
        push!(hyperedge_dist[panel][s], mean(final_hyperedges[s]))
        push!(hyperedge_dist_error[panel][s], std(final_hyperedges[s]))
    end
end

for panel in panels
    for s in 2:max_size
        scatter!(axes_matrix[panel], p_sweep, hyperedge_dist[panel][s];
                 markersize=8, marker=:xcross, label="simulation",
                 color=linecolors[s - 1])
        if draw_error_bars
            errorbars!(axes_matrix[panel], p_sweep, hyperedge_dist[panel][s],
                       hyperedge_dist_error[panel][s];
                       whiskerwidth=7, linewidth=0.5, color=linecolors[s - 1])
        end
    end
end

for panel in panels
    # compute the analytical solution
    params = load_params(joinpath(input_folder, "batch_001", "input_params.json"))
    tspan = (0.0, 100.0)

    rules = split("$panel", "_")
    params.model_params.propagation_rule = rules[1] == "maj" ? MajorityVoting() :
                                           ProportionalVoting()
    params.model_params.adaptivity_rule = rules[2] == "rts" ? RewireToSame() :
                                          RewireToRandom()

    function compute_hyperedge_dist(p)
        params.model_params.adaptivity_prob = p
        t, sol = moment_expansion(params, tspan, moment_closure)

        hyperedge_dist = zeros(max_size - 1)
        for motif in all_motifs(max_size)
            if AH.order(motif) != 1
                continue
            end
            size = AH.size(motif)
            hyperedge_dist[size - 1] += sol[motif][end]
        end

        return hyperedge_dist
    end

    p_analytical, hyperedge_dist = adaptive_sweep((0.0, 1.0),
                                                  compute_hyperedge_dist,
                                                  max_Δp, max_Δh)

    p_analytical = Vector{Float64}(p_analytical)
    # transpose the vector of vector of solutions
    hyperedge_dist = hcat(hyperedge_dist...)

    # draw the analytical solution
    for size in 2:max_size
        lines!(axes_matrix[panel], p_analytical, hyperedge_dist[size - 1, 1:end];
               color=linecolors[size - 1],
               linewidth=1.7, linestyle=:dash, label="mean-field")
    end
end

# draw the legend
plot_types = [LineElement(; color=linecolors[1], linewidth=1.7, linestyle=:dash),
              MarkerElement(; color=linecolors[1], marker=:xcross)]
colors = [PolyElement(; color=linecolors[i]) for i in 1:(max_size - 1)]
plot_type_labels = ["mean-field", "simulation"]
color_labels = ["size $size" for size in 2:max_size]

leg = Legend(fig[1, 2],
             [plot_types, colors],
             [plot_type_labels, color_labels],
             ["Plot type:", "Color:"])
leg.gridshalign = :left
leg.gridsvalign = :top

display(fig)

if prompt
    prompt_for_save(filename, fig)
end