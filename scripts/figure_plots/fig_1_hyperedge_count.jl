using AdaptiveHypergraphs
const AH = AdaptiveHypergraphs

include("figure_plotting_tools.jl")

# ======================================================
# -------------------- INPUT ---------------------------

input_folder = joinpath(projectdir(),
                        "./data/run_2023-04-14_15-08-19_motifs_all_rules")

title = "Total number of hyperedges"

num_nodes = 10000
max_size = 4
hyperedge_colormap = :thermal
linecolors = get(colorschemes[hyperedge_colormap], [1, 2.5, 3], (1, max_size))

# show a prompt whether the figure should be saved
prompt = true

filename = "./figures/fig_1_hyperedge_count.pdf"

batch_labels = ["Prop., rewire-to-random", "Majority, rewire-to-random",
                "Prop., rewire-to-same",
                "Majority, rewire-to-same"]

# ------------------------------------------------------
# ======================================================

# create figure
size_inches = (3.40, 2.8)
size_pt = 1.5 * 72 .* size_inches
fig = Figure(; resolution=size_pt, figure_padding=(7.0, 10.0, 5.0, 5.0))
subplots = fig[1, 1] = GridLayout()

axes_matrix = []

for row in 1:2
    push!(axes_matrix, [])
    for col in 1:2
        ax = Axis(subplots[row, col];
                  titlesize=12,
                  xticksize=3,
                  yticksize=3,
                  ylabelpadding=2,
                  xticklabelsize=10,
                  yticklabelsize=10,
                  xticklabelspace=10.0,
                  xticklabelsvisible=true,
                  yticklabelsvisible=true,
                  xlabelsize=12,
                  ylabelsize=12,
                  titlefont=:regular)
        ax.xticks = 0:20:100
        ax.yticks = 0:2500:12000
        if row == 1
            ax.xticksvisible = false
            ax.xticklabelsvisible = false
        else
            ax.xlabel = L"$t$"
        end
        if col == 2
            ax.yticksvisible = false
            ax.yticklabelsvisible = false
        else
            ax.ylabel = L"$\mathbf{M}(t)$"
        end

        push!(axes_matrix[row], ax)
    end
end

colgap!(subplots, 15)
rowgap!(subplots, 15)

data_folder = DataFolder(input_folder, :simulation)

for (batch_folder, batch_num) in data_folder
    row = (batch_num + 1) ÷ 2
    col = mod1(batch_num, 2)

    axes_matrix[row][col].title = batch_labels[batch_num]

    num_hyperedges_max = Dict(size => Float64[] for size in 2:max_size)
    num_hyperedges_min = Dict(size => Float64[] for size in 2:max_size)
    time = Float64[]

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

    for size in 2:max_size
        band!(axes_matrix[row][col], time, num_hyperedges_min[size],
              num_hyperedges_max[size];
              color=(linecolors[size - 1], 0.3), label="simulation")
    end

    xlims!(axes_matrix[row][col], (0, 100))
    ylims!(axes_matrix[row][col], (0, 12400))

    # find the analytical solution
    params = load_params(joinpath(input_folder, "batch_003", "input_params.json"))
    tspan = (0.0, 100.0)
    t, sol = moment_expansion(params, tspan, moment_closure)
    total_hyperedges_analytical = Dict(s => zero(t) for s in 2:max_size)
    for motif in all_motifs(max_size)
        if AH.order(motif) != 1
            continue
        end
        size = AH.size(motif)
        total_hyperedges_analytical[size] .+= sol[motif]
    end

    for size in 2:max_size
        lines!(axes_matrix[row][col], t, total_hyperedges_analytical[size];
               color=linecolors[size - 1],
               linewidth=1.7, linestyle=:dash, label="mean-field, size $size")
    end
end

plot_types = [LineElement(; color=linecolors[1], linewidth=1.7, linestyle=:dash),
              PolyElement(; color=(linecolors[1], 0.3))]
colors = [PolyElement(; color=linecolors[i]) for i in 1:(max_size - 1)]
plot_type_labels = ["mean-field", "simulation"]
color_labels = ["size $size" for size in 2:max_size]

"""
leg = Legend(fig[2, 1],
             [plot_types, colors],
             [plot_type_labels, color_labels],
             ["Line type:", "Color:"];
             orientation=:horizontal,
             labelsize=12,
             titlesize=12,
             titlefont="Latin Modern Roman Bold")
leg.gridshalign = :left
leg.gridsvalign = :top
leg.nbanks = 3
"""

display(fig)

if prompt
    prompt_for_save(filename, fig; pt_per_unit=0.666)
end