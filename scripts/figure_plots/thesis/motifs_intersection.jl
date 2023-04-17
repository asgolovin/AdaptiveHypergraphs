using StatsBase

include("figure_plotting_tools.jl")

#using GLMakie
#GLMakie.activate!()

using AdaptiveHypergraphs
const AH = AdaptiveHypergraphs

# ======================================================
# -------------------- INPUT ---------------------------

input_folder = joinpath(projectdir(),
                        "data/thesis/run_2023-03-12_16-02-12_motif_intersection_N_1e3")
data_folder = DataFolder(input_folder, :simulation)

num_nodes = 1000
max_size = 4
motifs = filter(x -> AH.order(x) == 2, all_motifs(max_size))
int_size_to_motifs = Dict(s => AbstractMotif[] for s in 1:max_size)
for s in 1:max_size
    filtered_motifs = filter(x -> x.int.A + x.int.B == s, motifs)
    append!(int_size_to_motifs[s], filtered_motifs)
end
node_colormap = :RdYlGn_6
hyperedge_colormap = :thermal

linecolors = get(colorschemes[hyperedge_colormap], [1, 2, 3, 4], (1, max_size + 1))

prompt = true

filename = "./figures/motifs_intersection.png"
ylabel = "# of motifs"
ylims = (0, 7000)
yticks = vcat([0], [10^i for i in 0:4])

# ------------------------------------------------------
# ======================================================

# create figure
size_inches = (6, 4)
size_pt = 1.5 * 72 .* size_inches
fig = Figure(; resolution=size_pt)

ax = Axis(fig[1, 1];
          yscale=Makie.pseudolog10,
          titlesize=15,
          title="Number of motifs with different intersection sizes",
          xticklabelsize=12,
          yticklabelsize=12,
          xlabelsize=13,
          ylabelsize=13,
          titlefont="Latin Modern Roman")

ax.xlabel = L"time $t$"
ax.ylabel = ylabel
ax.yticks = yticks

for (batch_folder, batch_num) in data_folder
    for (run_folder, run_num) in batch_folder
        if run_num != 1
            continue
        end
        indices, values = load_data("state_count_A.csv", run_folder)
        int_size_max = Dict(s => zero(indices) for s in 1:max_size)
        int_size_min = Dict(s => zero(indices) .+ 1e9 for s in 1:max_size)
        for int_size in 1:max_size
            for motif in int_size_to_motifs[int_size]
                indices, values = load_data("motif_count_$motif.csv", run_folder)
                int_size_max[int_size] .= max.(int_size_max[int_size], values)
                int_size_min[int_size] .= min.(int_size_min[int_size], values)
            end
        end
        for int_size in 1:max_size
            band!(ax, indices, int_size_min[int_size], int_size_max[int_size];
                  color=(linecolors[int_size], 0.3),
                  label="intersection size = $int_size")
            if int_size != 1
                lines!(ax, indices, int_size_max[int_size] ./ int_size_max[1];
                       color=linecolors[int_size])
            end
        end
    end
end

Legend(fig[1, 2], ax)

display(fig)

if prompt
    prompt_for_save(filename, fig; pt_per_unit=0.666)
end