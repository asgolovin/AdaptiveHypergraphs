using StatsBase

plot_analytical = true

if plot_analytical
    using AdaptiveHypergraphs
    const AH = AdaptiveHypergraphs
end

include("figure_plotting_tools.jl")

# ======================================================
# -------------------- INPUT ---------------------------

input_folder = joinpath(projectdir(),
                        "final_data/run_2023-02-22_16-44-20_motifs_D3")
data_folder = DataFolder(input_folder, :simulation)

num_nodes = 10000
max_size = 3
labels = filter(x -> AH.order(x) == 1, all_motifs(max_size))
active_labels = filter(x -> x.left[AH.A] != 0 && x.left[AH.B] != 0, labels)
inactive_labels = filter(x -> x ∉ active_labels, labels)
node_colormap = :RdYlGn_6
hyperedge_colormap = :thermal

# colorscheme = colorschemes[node_colormap]
# linecolors = Dict()
# linestyles = Dict()
# for (i, label) in enumerate(labels)
#     AB_ratio = label.left_total[AH.A] / (size(label)[1])
#     linecolors[label] = get(colorscheme, AB_ratio, (0.0, 1.0))
#     linestyles[label] = [nothing, :dash, :dot][size(label)[1]]
# end

linecolors = get(colorschemes[hyperedge_colormap], [1, 2.5, 3], (1, max_size + 1))

prompt = true

filename = "./figures/motifs_analytical.pdf"
ylabel = "# of motifs"
ylims = (0, 7000)
yticks = 0:2000:8000

# ------------------------------------------------------
# ======================================================

# create figure
size_inches = (6, 4)
size_pt = 1.5 * 72 .* size_inches
fig = Figure(; resolution=size_pt)

ax_inactive = Axis(fig[1, 1];
                   titlesize=15,
                   title="Inactive",
                   xticklabelsize=12,
                   yticklabelsize=12,
                   xlabelsize=13,
                   ylabelsize=13,
                   titlefont="Latin Modern Roman")
ax_active = Axis(fig[2, 1];
                 titlesize=15,
                 title="Active",
                 xticklabelsize=12,
                 yticklabelsize=12,
                 xlabelsize=13,
                 ylabelsize=13,
                 titlefont="Latin Modern Roman")
hidexdecorations!(ax_inactive; grid=false)
ax_active.xlabel = L"time $t$"
ax_inactive.ylabel = ylabel
ax_active.ylabel = ylabel
ax_active.yticks = yticks
ax_inactive.yticks = yticks
xlims!(ax_inactive, (0, 50))
xlims!(ax_active, (0, 50))
ylims!(ax_inactive, ylims)
ylims!(ax_active, ylims)

indices, values = load_data("state_count_A.csv",
                            joinpath(input_folder, "batch_001", "run_001"))
motif_max = Dict(label => zero(indices) for label in labels)
motif_min = Dict(label => zero(indices) .+ 1e9 for label in labels)

for (batch_folder, batch_num) in data_folder
    for (run_folder, run_num) in batch_folder
        for label in labels
            indices, values = load_data("motif_count_$label.csv", run_folder)

            motif_max[label] .= max.(motif_max[label], values)
            motif_min[label] .= min.(motif_min[label], values)
        end
    end
end

motif_max[Label("[A2]")] .= max.(motif_max[Label("[A2]")], motif_max[Label("[B2]")])
motif_min[Label("[A2]")] .= min.(motif_min[Label("[A2]")], motif_min[Label("[B2]")])
motif_max[Label("[A3]")] .= max.(motif_max[Label("[A3]")], motif_max[Label("[B3]")])
motif_min[Label("[A3]")] .= min.(motif_min[Label("[A3]")], motif_min[Label("[B3]")])
motif_max[Label("[A2B]")] .= max.(motif_max[Label("[A2B]")], motif_max[Label("[AB2]")])
motif_min[Label("[A2B]")] .= min.(motif_min[Label("[A2B]")], motif_min[Label("[AB2]")])

deleteat!(labels, findfirst(x -> x == Label("[B2]"), labels))
deleteat!(labels, findfirst(x -> x == Label("[B3]"), labels))
deleteat!(labels, findfirst(x -> x == Label("[AB2]"), labels))

# plot the analytical results
params = load_params(joinpath(input_folder, "batch_001", "input_params.json"))
tspan = (0.0, 50.0)
t, sol = moment_expansion(params, tspan, moment_closure)

label_labels = Dict(AH.Label("[A2]") => L"$[A^2]$ and $[B^2]$",
                    AH.Label("[A3]") => L"$[A^3]$ and $[B^3]$",
                    AH.Label("[AB]") => L"$[AB]$",
                    AH.Label("[A2B]") => L"$[A^2B]$ and $[AB^2]$")

for label in labels
    size = AH.size(label)[1]
    ax = label ∈ active_labels ? ax_active : ax_inactive
    if label.left[AH.A] >= label.left[AH.B]
        lines!(ax, t, sol[label];
               linewidth=1.7,
               linestyle=:dash,
               color=linecolors[size - 1], label=label_labels[label])
    end
end

for label in labels
    ax = label ∈ active_labels ? ax_active : ax_inactive
    size = AH.size(label)[1]
    println("$label, $size")
    if size == 2
        band!(ax, indices, motif_min[label], motif_max[label];
              color=(linecolors[size - 1], 0.3),
              label="simulation")
    else
        band!(ax, indices, motif_min[label], motif_max[label];
              color=(linecolors[size - 1], 0.3))
    end
end

Legend(fig[1, 2], ax_inactive; merge=true, orientation=:vertical, labelsize=12)
Legend(fig[2, 2], ax_active; merge=true, orientation=:vertical, labelsize=12)

display(fig)

if prompt
    prompt_for_save(filename, fig; pt_per_unit=0.666)
end