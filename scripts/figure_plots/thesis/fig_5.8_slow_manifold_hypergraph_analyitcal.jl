using AdaptiveHypergraphs
const AH = AdaptiveHypergraphs

include("figure_plotting_tools.jl")

# ======================================================
# -------------------- INPUT ---------------------------

#input_folder = joinpath(projectdir(),
#"data/thesis/#run_2023-02-06_21-50-23_prop_voting_rtr")
input_folder = joinpath(projectdir(),
                        "data/thesis/run_2023-02-23_17-46-21_slow_manifold_p08")

title = "Slow manifold plot"

num_nodes = 10000
max_size = 4
hyperedge_colormap = :thermal
linecolors = get(colorschemes[hyperedge_colormap], [1, 2.5, 3], (1, max_size))

prompt = true

filename = "./figures/thesis/slow_manifold_hypergraph_analytical.pdf"

# ------------------------------------------------------
# ======================================================

# create figure
size_inches = (5, 3)
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
ax.xlabel = L"magnetization $m$"
ax.ylabel = L"density of active links $\rho_{d_i}$"
#ax.yticks = 0.0:0.1:0.6
ax.xticks = -1:0.5:1

skip_initial = 200

magnetization_data = []
density_data = Dict(s => [] for s in 2:max_size)

data_folder = DataFolder(input_folder, :simulation)

indices, state_count_A = load_data("state_count_A.csv",
                                   joinpath(input_folder, "batch_001", "run_001"))
magnetization_bins = collect(-1:0.01:1)
hyperedge_density_max = Dict(s => zero(magnetization_bins) for s in 2:max_size)
hyperedge_density_min = Dict(s => zero(magnetization_bins) .+ 1e9 for s in 2:max_size)

for (batch_folder, batch_num) in data_folder
    println("$batch_num")
    for (run_folder, run_num) in batch_folder
        indices, state_count_A = load_data("state_count_A.csv", run_folder)
        magnetization = 2 .* (state_count_A ./ num_nodes) .- 1
        magnetization = magnetization[skip_initial:end]
        append!(magnetization_data, magnetization)

        # the index of the bin for every value in magnetization
        bins = Int64.(round.(magnetization .* 100)) .+ 101

        for size in 2:max_size
            filename = "active_hyperedge_count_$size.csv"
            indices, active_hyperedeges = load_data(filename, run_folder)
            active_hyperedeges = active_hyperedeges[skip_initial:end]
            active_density = active_hyperedeges ./ num_nodes
            append!(density_data[size], active_density)

            for i in 1:length(magnetization_bins)
                max_density = maximum(active_density[bins .== i]; init=0)
                min_density = minimum(active_density[bins .== i]; init=1e9)
                hyperedge_density_max[size][i] = max(hyperedge_density_max[size][i],
                                                     max_density)
                hyperedge_density_min[size][i] = min(hyperedge_density_min[size][i],
                                                     min_density)
            end
        end
    end
end

for size in 2:max_size
    hyperedge_density_min[size][hyperedge_density_min[size] .> 1e8] .= 0.0

    if size == 3
        band!(ax, magnetization_bins, hyperedge_density_min[size],
              hyperedge_density_max[size];
              color=(linecolors[size - 1], 0.3),
              label="simulation")
    else
        band!(ax, magnetization_bins, hyperedge_density_min[size],
              hyperedge_density_max[size];
              color=(linecolors[size - 1], 0.3))
    end
end

# x = magnetization_bins
# for size in 2:max_size
#     p, cov = fit_parabola(magnetization_data,
#                           density_data[size], false)
#     lines!(ax, x, p.(x); color=linecolors[size - 1], label="polynomial fit")
# end

xlims!(ax, (-1, 1))
ylims!(ax, (0, 0.3))

magn_analytical = collect(-1:0.1:1)
density_analytical = Dict(s => zero(magn_analytical) for s in 2:max_size)

params = load_params(joinpath(input_folder, "batch_001", "input_params.json"))
params.model_params.adaptivity_prob = 0.85

for (i, m0) in enumerate(magn_analytical)
    println("$m0")
    params.network_params.state_A_prob = (m0 + 1) * 0.5
    tspan = (0.0, 1000.0)
    t, sol = moment_expansion(params, tspan, moment_closure)
    for label in all_motifs(max_size)
        if AH.order(label) == 1 && label.left[AH.A] > 0 && label.left[AH.B] > 0
            size = AH.size(label)[1]
            density_analytical[size][i] += sol[label][end] / num_nodes
        end
    end
end

for size in 2:max_size
    lines!(ax, magn_analytical, density_analytical[size];
           color=linecolors[size - 1],
           linestyle=:dash, linewidth=1.7, label="mean-field, size $size")
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