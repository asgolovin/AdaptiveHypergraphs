using StatsBase

using AdaptiveHypergraphs
const AH = AdaptiveHypergraphs

include("figure_plotting_tools.jl")

# ======================================================
# -------------------- INPUT ---------------------------

draw_error_bars = false

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
              :prop_rts => "Proportional voting,\nrewire-to-same",
              :prop_rtr => "Proportional voting,\nrewire-to-random")

input_folders = Dict(:maj_rts => folder_maj_rts,
                     :maj_rtr => folder_maj_rtr,
                     :prop_rts => folder_prop_rts,
                     :prop_rtr => folder_prop_rtr)

num_nodes = 10000
max_size = 4
hyperedge_colormap = :thermal
linecolors = get(colorschemes[hyperedge_colormap], [1, 2.5, 3], (1, max_size))

prompt = false

# remove first seconds of the simulation
initial_cutoff = 30

filename = "./figures/fig_4_p_sweep.pdf"

xlabel = L"$p$"
ylabel = L"$\rho_i$"
xticks = 0:0.2:1
yticks = 0:0.2:1
xlims = (-0.02, 1.02)
ylims = (nothing, 0.5)

max_Δp = 1
max_Δh = 0.15

cache_file = "./analytical_stable_points.csv"

# ------------------------------------------------------
# ======================================================

function sim_parabola_height(batch_folder)
    # Collect data from all runs, then fit a parabola to all the data
    magnetization_all_runs = Float64[]
    active_hyperedges_all_runs = Dict(s => Float64[] for s in 2:max_size)

    params = load_params(joinpath(batch_folder.folder, "input_params.json"))
    p = params.model_params.adaptivity_prob

    for (run_folder, _) in batch_folder
        indices, state_count_A = load_data("state_count_A.csv", run_folder)
        mag = 2 .* state_count_A ./ num_nodes .- 1

        mask = indices .> initial_cutoff
        # filter out all values which lie in the center
        roots = abs(mag[end])
        mask = mask .| (abs.(mag) .> roots * 0.8)

        # if there are too few data points, take only the last value
        if sum(mask) < 10
            mask = fill(false, length(mag))
            mask[end] = true
        end
        append!(magnetization_all_runs, mag[mask])

        for size in 2:max_size
            _, active_hyperedges = load_data("active_hyperedge_count_$size.csv", run_folder)
            active_hyperedges = Vector{Float64}(active_hyperedges)
            active_hyperedges ./= num_nodes
            append!(active_hyperedges_all_runs[size], active_hyperedges[mask])
        end
    end

    height_by_size = Dict(s => 0.0 for s in 2:max_size)
    height_error = Dict(s => 0.0 for s in 2:max_size)

    for size in 2:max_size
        parabola, cov = fit_parabola(magnetization_all_runs,
                                     active_hyperedges_all_runs[size])
        height_by_size[size] = parabola(0.0)
        height_error[size] = sqrt(cov[1, 1])
    end
    return p, height_by_size, height_error
end

function analytical_parabola_height(params)
    tspan = (0.0, 200.0)
    params.network_params.state_A_prob = 0.5
    num_nodes = params.network_params.num_nodes

    function compute_parabola_height(p)
        params.model_params.adaptivity_prob = p
        params.network_params.state_A_prob = 0.5
        t, sol = moment_expansion(params, tspan, moment_closure)

        height_by_size = Dict(s => 0.0 for s in 2:max_size)
        for a in 1:(max_size - 1), b in 1:(max_size - a)
            active_hyperedges = sol[OrderOneMotif(a, b)][end]
            active_hyperedges /= num_nodes
            height_by_size[a + b] += active_hyperedges
        end

        return height_by_size
    end

    p_sweep = vcat([0.0, 0.01, 0.02], collect(0.1:0.1:1.0))
    height_sweep = []
    for p in p_sweep
        height = compute_parabola_height(p)
        push!(height_sweep, height)
    end

    return p_sweep, height_sweep
end

# create figure
fig = create_figure(:large, 3 / 2)
ax_gridpos = fig[1, 1]

axis = add_four_rules_axes!(fig, titles, xlabel, ylabel, xticks, yticks, xlims, ylims)

for panel in panels
    p_sweep = Float64[]
    parabola_height = Dict(s => Float64[] for s in 2:max_size)
    parabola_height_error = Dict(s => Float64[] for s in 2:max_size)

    data_folder = DataFolder(input_folders[panel], :simulation)
    for (batch_folder, batch_num) in data_folder
        p, height, error = sim_parabola_height(batch_folder)
        push!(p_sweep, p)
        for size in 2:max_size
            push!(parabola_height[size], height[size])
            push!(parabola_height_error[size], error[size])
        end
    end

    for size in 2:max_size
        scatter!(axis[panel], p_sweep, parabola_height[size];
                 markersize=8, marker=:xcross, label="simulation",
                 color=linecolors[size - 1])
        if draw_error_bars
            errorbars!(axis[panel], p_sweep, parabola_height[size],
                       parabola_height_error[size];
                       whiskerwidth=7, linewidth=0.5)
        end
    end

    params = load_params(joinpath(input_folders[panel], "batch_001", "input_params.json"))
    p, parabola_height = analytical_parabola_height(params)
    for size in 2:max_size
        height = [d[size] for d in parabola_height]
        lines!(axis[panel], p, height;
               linewidth=1.7,
               linestyle=:dash, label="mean-field",
               color=linecolors[size - 1])
    end
end

plot_types = [LineElement(; color=linecolors[1], linewidth=1.7, linestyle=:dash),
              MarkerElement(; marker=:xcross, markersize=8, color=linecolors[1])]
colors = [PolyElement(; color=linecolors[i]) for i in 1:(max_size - 1)]
plot_type_labels = ["mean-field", "simulation"]
color_labels = ["size $size" for size in 2:max_size]

leg = Legend(fig[1, 2],
             [plot_types, colors],
             [plot_type_labels, color_labels],
             ["Line type", "Color:"];
             orientation=:vertical,
             labelsize=12,
             titlesize=12)
leg.gridshalign = :left

display(fig)

if prompt
    prompt_for_save(filename, fig)
end