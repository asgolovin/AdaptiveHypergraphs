using AdaptiveHypergraphs
const AH = AdaptiveHypergraphs

include("figure_plotting_tools.jl")

# ======================================================
# -------------------- INPUT ---------------------------

input_folder = joinpath(projectdir(),
                        "results/run_2023-06-16_10-40-16_motifs_maj")

num_nodes = 10000
max_size = 4
linecolors = hyperedge_linecolors(max_size)

# show a prompt whether the figure should be saved
prompt = false

filename = "./figures/fig_1_hyperedge_count.pdf"

panels = [:maj_rts, :maj_rtr, :prop_rts, :prop_rtr]

titles = Dict(:maj_rts => "Majority voting, rewire-to-same",
              :maj_rtr => "Majority voting, rewire-to-random",
              :prop_rts => "Prop. voting, rewire-to-same",
              :prop_rtr => "Prop. voting, rewire-to-random")

batch_to_panel = Dict(1 => :prop_rtr,
                      2 => :maj_rtr,
                      3 => :prop_rts,
                      4 => :maj_rts)

xlabel = L"$t$"
ylabel = L"$\mathbf{M}(t)$"
xticks = 0:20:100
yticks = 0:2500:12000
xlims = (0.0, 100.0)
ylims = (0, 12400)

# ------------------------------------------------------
# ======================================================

# create figure
fig = create_figure(:large, 2 / 1)
subplots = fig[1, 1] = GridLayout()

axes_matrix = add_four_rules_axes!(subplots, titles,
                                   xlabel, ylabel,
                                   xticks, yticks,
                                   xlims, ylims)

data_folder = DataFolder(input_folder, :simulation)

for (batch_folder, batch_num) in data_folder
    panel = batch_to_panel[batch_num]

    num_hyperedges_max = Dict(size => Float64[] for size in 2:max_size)
    num_hyperedges_min = Dict(size => Float64[] for size in 2:max_size)
    time = Float64[]

    # accumulate maximum and minimum number of hyperedges
    for (run_folder, run_num) in batch_folder
        for size in 2:max_size
            indices, values = load_data("hyperedge_count_$size.csv", run_folder)
            num_hyperedges_max[size] = safe_max(num_hyperedges_max[size], values)
            num_hyperedges_min[size] = safe_min(num_hyperedges_min[size], values)
            if length(indices) > length(time)
                time = vcat(time, indices[(length(time) + 1):end])
            end
        end
    end

    # draw the simulation results
    for size in 2:max_size
        band!(axes_matrix[panel], time, num_hyperedges_min[size],
              num_hyperedges_max[size];
              color=(linecolors[size - 1], 0.5))
    end

    # find the analytical solution
    params = load_params(joinpath(batch_folder.folder, "input_params.json"))
    tspan = xlims
    t, sol = moment_expansion(params, tspan, moment_closure)
    total_hyperedges_analytical = Dict(s => zero(t) for s in 2:max_size)
    for motif in all_motifs(max_size)
        if AH.order(motif) != 1
            continue
        end
        size = AH.size(motif)
        total_hyperedges_analytical[size] .+= sol[motif]
    end

    # draw the analytical solution
    for size in 2:max_size
        lines!(axes_matrix[panel], t, total_hyperedges_analytical[size];
               color=linecolors[size - 1],
               linewidth=1.7, linestyle=:dash)
    end
end

# draw the legend
plot_types = [LineElement(; color=linecolors[1], linewidth=1.7, linestyle=:dash),
              PolyElement(; color=(linecolors[1], 0.3))]
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