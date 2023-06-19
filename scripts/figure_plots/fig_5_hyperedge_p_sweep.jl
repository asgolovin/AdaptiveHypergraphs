using StatsBase

using AdaptiveHypergraphs
const AH = AdaptiveHypergraphs

include("figure_plotting_tools.jl")

# ======================================================
# -------------------- INPUT ---------------------------

draw_error_bars = true

folder_maj_rts = joinpath(projectdir(),
                          "data/run_2023-04-24_14-46-47_maj_voting_rts")
folder_maj_rtr = joinpath(projectdir(),
                          "data/run_2023-04-24_14-47-49_maj_voting_rtr")
folder_prop_rts = joinpath(projectdir(),
                           "data/run_2023-05-11_14-55-38_prop_voting_rts")
folder_prop_rtr = joinpath(projectdir(),
                           "data/run_2023-04-26_06-09-49_prop_voting_rtr")

panels = [:maj_rts, :maj_rtr, :prop_rts, :prop_rtr]
titles = Dict(:maj_rts => "Majority voting,\nrewire-to-same",
              :maj_rtr => "Majority voting,\nrewire-to-random",
              :prop_rts => "Proportional voting\nrewire-to-same",
              :prop_rtr => "Proportional voting\nrewire-to-random")

input_folders = Dict(:maj_rts => folder_maj_rts,
                     :maj_rtr => folder_maj_rtr,
                     :prop_rts => folder_prop_rts,
                     :prop_rtr => folder_prop_rtr)

num_nodes = 10000
max_size = 4
linecolors = hyperedge_linecolors(max_size)

prompt = false

# number of active hyperedges below which the simulation counts as converged
converged_threshold = 10

filename = "./figures/fig_5_hyperedges_p_sweep.pdf"

xlabel = L"$p$"
ylabel = L"\mathbf{M}(t)"
xticks = 0:0.2:1
yticks = 0:2500:12000
xlims = (-0.02, 1.02)
ylims = (0, 12500)

max_Δp = 0.5
max_Δh = 800

# ------------------------------------------------------
# ======================================================

# create figure
fig = create_figure(:large, 3 / 2)
ax_gridpos = fig[1, 1] = GridLayout()

axes_matrix = add_four_rules_axes!(ax_gridpos, titles,
                                   xlabel, ylabel,
                                   xticks, yticks,
                                   xlims, ylims)

for panel in panels
    p_sweep = Float64[]
    hyperedge_dist = Dict(s => Float64[] for s in 2:max_size)
    hyperedge_dist_error = Dict(s => Float64[] for s in 2:max_size)

    data_folder = DataFolder(input_folders[panel], :simulation)
    for (batch_folder, batch_num) in data_folder
        params = load_params(joinpath(batch_folder.folder, "input_params.json"))
        p = params.model_params.adaptivity_prob
        push!(p_sweep, p)

        final_hyperedges = Dict(s => Float64[] for s in 2:max_size)
        for (run_folder, run_num) in batch_folder
            for s in 2:max_size
                time, hyperedge_count = load_data("hyperedge_count_$s.csv", run_folder)

                # take all solutions after 100 seconds (or the last 10% if the simulation lasted less than 100 seconds) and average them
                if time[end] < 100
                    mask = time .> time[end] * 0.9
                else
                    mask = time .> 100
                end
                push!(final_hyperedges[s], mean(hyperedge_count[mask]))
            end
        end

        for s in 2:max_size
            push!(hyperedge_dist[s], mean(final_hyperedges[s]))
            push!(hyperedge_dist_error[s], std(final_hyperedges[s]))
        end
    end

    for s in 2:max_size
        scatter!(axes_matrix[panel], p_sweep, hyperedge_dist[s];
                 markersize=8, marker=:xcross, label="simulation",
                 color=linecolors[s - 1])
        if draw_error_bars
            errorbars!(axes_matrix[panel], p_sweep, hyperedge_dist[s],
                       hyperedge_dist_error[s];
                       whiskerwidth=7, linewidth=0.5, color=linecolors[s - 1])
        end
    end

    # compute the analytical solution
    params = load_params(joinpath(input_folders[panel], "batch_001", "input_params.json"))
    tspan = (0.0, 100.0)

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