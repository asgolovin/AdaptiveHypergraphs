using AdaptiveHypergraphs
const AH = AdaptiveHypergraphs

include("figure_plotting_tools.jl")

# ======================================================
# -------------------- INPUT ---------------------------

input_folder = joinpath(projectdir(),
                        "data/thesis/run_2023-02-06_21-50-23_prop_voting_rtr")

title = "Total number of hyperedges"

num_nodes = 10000
max_size = 4
hyperedge_colormap = :thermal
linecolors = get(colorschemes[hyperedge_colormap], [1, 2.5, 3], (1, max_size))

# plot every nth point 
skip_middle = 10

prompt = true

# which run to plot in bold
highlight_run = 5

filename = "./figures/thesis/hyperedge_dist_analytical.pdf"

# ------------------------------------------------------
# ======================================================

# create figure
size_inches = (4.5, 2.8)
size_pt = 1.5 * 72 .* size_inches
fig = Figure(; resolution=size_pt)

ax = Axis(fig[1, 1];
          titlesize=15,
          title=title,
          xticklabelsize=12,
          yticklabelsize=12,
          xlabelsize=13,
          ylabelsize=13,
          titlefont="Latin Modern Roman")
ax.xticks = 0:20:100
ax.yticks = 0:2500:12000

ax.xlabel = L"time $t$"
ax.ylabel = "# of hyperedges"

data_folder = DataFolder(input_folder, :simulation)

indices, values = load_data("state_count_A.csv",
                            joinpath(input_folder, "batch_003", "run_001"))
num_hyperedges_max = Dict(size => zero(indices) for size in 2:max_size)
num_hyperedges_min = Dict(size => zero(indices) .+ 1e9 for size in 2:max_size)

for (batch_folder, batch_num) in data_folder
    if batch_num != 3
        continue
    end
    for (run_folder, run_num) in batch_folder
        for size in 2:max_size
            indices, values = load_data("hyperedge_count_$size.csv", run_folder)

            num_hyperedges_max[size] .= max.(num_hyperedges_max[size], values)
            num_hyperedges_min[size] .= min.(num_hyperedges_min[size], values)
        end
    end
end

for size in 2:max_size
    if size == 3
        band!(ax, indices, num_hyperedges_min[size], num_hyperedges_max[size];
              color=(linecolors[size - 1], 0.3), label="simulation")
    else
        band!(ax, indices, num_hyperedges_min[size], num_hyperedges_max[size];
              color=(linecolors[size - 1], 0.3))
    end
end

xlims!(ax, (0, 100))
ylims!(ax, (0, 12400))

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
    lines!(ax, t, total_hyperedges_analytical[size]; color=linecolors[size - 1],
           linewidth=1.7, linestyle=:dash, label="mean-field, size $size")
end

plot_types = [LineElement(; color=linecolors[1], linewidth=1.7, linestyle=:dash),
              PolyElement(; color=(linecolors[1], 0.3))]
colors = [PolyElement(; color=linecolors[i]) for i in 1:(max_size - 1)]
plot_type_labels = ["mean-field", "simulation"]
color_labels = ["size $size" for size in 2:max_size]

leg = Legend(fig[1, 2],
             [plot_types, colors],
             [plot_type_labels, color_labels],
             ["", "Color:"];
             orientation=:vertical,
             labelsize=12,
             titlesize=12,
             titlefont="Latin Modern Roman Bold")
leg.gridshalign = :left

display(fig)

if prompt
    prompt_for_save(filename, fig; pt_per_unit=0.666)
end