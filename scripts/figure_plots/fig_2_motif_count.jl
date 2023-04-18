using AdaptiveHypergraphs
const AH = AdaptiveHypergraphs

include("figure_plotting_tools.jl")

# ======================================================
# -------------------- INPUT ---------------------------

input_folder = joinpath(projectdir(),
                        "./data/run_2023-04-14_15-08-19_motifs_all_rules")

title = "Total number of motifs"

motifs = filter(x -> typeof(x) <: OrderOneMotif, all_motifs(max_size))

num_nodes = 10000
max_size = 4
hyperedge_colormap = :thermal
node_colormap = :RdYlGn_6
#linecolors = get(colorschemes[hyperedge_colormap], [1, 2.5, 3], (1, max_size))
colorscheme = colorschemes[node_colormap]
linecolors = Dict()
for (i, motif) in enumerate(motifs)
    AB_ratio = motif.A / size(motif)
    linecolors[motif] = get(colorscheme, AB_ratio, (0.0, 1.0))
end

# show a prompt whether the figure should be saved
prompt = true

filename = "./figures/fig_2_motif_count.pdf"

batch_labels = ["Prop. voting,\n rewire-to-random", "Majority voting,\n rewire-to-random",
                "Prop. voting,\n rewire-to-same",
                "Majority voting,\n rewire-to-same"]

# every row has its own ymax and ticks
ymax = [11000, 5000, 1500]
yticks = [0:3000:11000, 0:2000:6000, 0:500:1500]

# ------------------------------------------------------
# ======================================================

# create figure
fig = create_figure(:large, 4 / 3)
subplots = fig[1, 1] = GridLayout()

axes_matrix = []

for row in 1:3
    push!(axes_matrix, [])
    for col in 1:4
        ax = Axis(subplots[row, col])
        ax.xticks = 0:20:100
        ax.yticks = yticks[row]
        if row != 3
            ax.xticksvisible = false
            ax.xticklabelsvisible = false
        else
            ax.xlabel = L"$t$"
        end
        if col != 1
            ax.yticksvisible = false
            ax.yticklabelsvisible = false
        else
            ax.ylabel = L"$\mathbf{M}(t)$"
        end

        push!(axes_matrix[row], ax)
    end
end

data_folder = DataFolder(input_folder, :simulation)

for (batch_folder, batch_num) in data_folder
    col = batch_num
    axes_matrix[1][col].title = batch_labels[batch_num]

    num_motifs_max = Dict(motif => Float64[] for motif in motifs)
    num_motifs_min = Dict(motif => Float64[] for motif in motifs)
    time = Float64[]

    for (run_folder, run_num) in batch_folder
        for motif in motifs
            indices, values = load_data("motif_count_$motif.csv", run_folder)
            num_motifs_max[motif] = safe_max(num_motifs_max[motif], values)
            num_motifs_min[motif] = safe_min(num_motifs_min[motif], values)
            if length(indices) > length(time)
                time = vcat(time, indices[(length(time) + 1):end])
            end
        end
    end

    for s in 2:max_size
        row = s - 1
        for motif in filter(x -> size(x) == s, motifs)
            band!(axes_matrix[row][col], time, num_motifs_min[motif],
                  num_motifs_max[motif];
                  color=(linecolors[motif], 0.5), label="simulation")
        end
        xlims!(axes_matrix[row][col], (0, 100))
        ylims!(axes_matrix[row][col], (-0.03 * ymax[row], ymax[row]))
    end

    # find the analytical solution
    params = load_params(joinpath(batch_folder.folder, "input_params.json"))
    tspan = (0.0, 100.0)
    t, sol = moment_expansion(params, tspan, moment_closure)

    for s in 2:max_size
        row = s - 1
        for motif in filter(x -> size(x) == s, motifs)
            lines!(axes_matrix[row][col], t, sol[motif];
                   color=linecolors[motif],
                   linewidth=1.7, linestyle=:dash, label="$motif")
        end
    end
end

# plot the legends for every row
for s in 2:max_size
    row = s - 1
    colors = [PolyElement(; color=linecolors[motif])
              for motif in filter(x -> size(x) == s, motifs)]
    color_labels = ["$motif" for motif in filter(x -> size(x) == s, motifs)]
    color_labels = [rpad(label, 10, " ") for label in color_labels]

    leg = Legend(subplots[row, 5],
                 colors,
                 color_labels;
                 orientation=:vertical,
                 labelsize=10,
                 titlesize=12,
                 width=80)
end

colgap!(subplots, 15)
rowgap!(subplots, 15)

display(fig)

if prompt
    prompt_for_save(filename, fig)
end