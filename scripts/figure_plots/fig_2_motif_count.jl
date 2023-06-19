using AdaptiveHypergraphs
const AH = AdaptiveHypergraphs

include("figure_plotting_tools.jl")

# ======================================================
# -------------------- INPUT ---------------------------

input_folder = joinpath(projectdir(),
                        "results/run_2023-06-16_10-40-16_motifs_maj")

num_nodes = 10000
max_size = 4
node_colormap = :RdYlGn_6

motifs = filter(x -> typeof(x) <: OrderOneMotif, all_motifs(max_size))

colorscheme = colorschemes[node_colormap]
linecolors = Dict()
for (i, motif) in enumerate(motifs)
    AB_ratio = motif.A / size(motif)
    linecolors[motif] = get(colorscheme, AB_ratio, (0.0, 1.0))
end

# show a prompt whether the figure should be saved
prompt = false

filename = "./figures/fig_2_motif_count.pdf"

panels = [:maj_rts, :maj_rtr, :prop_rts, :prop_rtr]

titles = Dict(:maj_rts => "Majority voting,\nrewire-to-same",
              :maj_rtr => "Majority voting,\nrewire-to-random",
              :prop_rts => "Prop. voting,\nrewire-to-same",
              :prop_rtr => "Prop. voting,\nrewire-to-random")

batch_to_panel = Dict(1 => :prop_rtr,
                      2 => :maj_rtr,
                      3 => :prop_rts,
                      4 => :maj_rts)

xlabel = L"$t$"
ylabel = "Number of motifs"
xticks = 0:20:100
xlims = (0.0, 100.0)
# every row has its own ymax and ticks
ymax = [14000, 5500, 2000]
yticks = [0:3000:14000, 0:2000:6000, 0:500:2000]

# ------------------------------------------------------
# ======================================================

# create figure
fig = create_figure(:large, 4 / 3)
subplots = fig[1, 1] = GridLayout()

axes_matrix = Dict(panel => [] for panel in panels)

for (col, panel) in enumerate(panels)
    for row in 1:(max_size - 1)
        local ax = Axis(subplots[row, col])
        ax.xticks = xticks
        ax.yticks = yticks[row]
        xlims!(ax, xlims)
        ylims!(ax, (-0.03 * ymax[row], ymax[row]))
        if row != 3
            ax.xticksvisible = false
            ax.xticklabelsvisible = false
        else
            ax.xlabel = xlabel
        end
        if col != 1
            ax.yticksvisible = false
            ax.yticklabelsvisible = false
        else
            ax.ylabel = ylabel
        end

        push!(axes_matrix[panel], ax)
    end
end

data_folder = DataFolder(input_folder, :simulation)

for (batch_folder, batch_num) in data_folder
    panel = batch_to_panel[batch_num]

    axes_matrix[panel][1].title = titles[panel]

    num_motifs_max = Dict(motif => Float64[] for motif in motifs)
    num_motifs_min = Dict(motif => Float64[] for motif in motifs)
    time = Float64[]

    # accumulate the maximum and minimum number of motifs
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

    # plot the motifs
    for s in 2:max_size
        row = s - 1
        for motif in filter(x -> size(x) == s, motifs)
            band!(axes_matrix[panel][row], time, num_motifs_min[motif],
                  num_motifs_max[motif];
                  color=(linecolors[motif], 0.5))
        end
    end

    # find the analytical solution
    params = load_params(joinpath(batch_folder.folder, "input_params.json"))
    tspan = xlims
    t, sol = moment_expansion(params, tspan, moment_closure)

    # plot the analytical solution
    for s in 2:max_size
        row = s - 1
        for motif in filter(x -> size(x) == s, motifs)
            lines!(axes_matrix[panel][row], t, sol[motif];
                   color=linecolors[motif],
                   linewidth=1.7, linestyle=:dash)
        end
    end
end

# plot the legends for every row
for s in 2:max_size
    row = s - 1
    local colors = [PolyElement(; color=linecolors[motif])
                    for motif in filter(x -> size(x) == s, motifs)]

    local color_labels = ["$motif" for motif in filter(x -> size(x) == s, motifs)]
    # replace the numbers by superscripts and add a tiny space ("hair space") after the superscripts
    color_labels = [replace(motif, "2" => "² ", "3" => "³ ", "4" => "⁴ ", " " => " ")
                    for motif in color_labels]
    # replace double hair space by a single hair space
    color_labels = [replace(motif, "  " => " ")
                    for motif in color_labels]

    local leg = Legend(subplots[row, 5],
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