using StatsBase

plot_analytical = true

if plot_analytical
    using AdaptiveHypergraphs
    const AH = AdaptiveHypergraphs
end

include("figure_plotting_tools.jl")

# ======================================================
# -------------------- INPUT ---------------------------

plot_type = :total # :active or :total
draw_error_bars = false

input_folder = joinpath(projectdir(),
                        "results/run_2023-02-08_03-12-17_p_sweep")
data_folder = DataFolder(input_folder, :simulation)

# time from which the trajectories converged to a stable distribution for each batch
starting_point = Dict(1 => 5,
                      2 => 50,
                      3 => 20,
                      4 => 20,
                      5 => 20,
                      6 => 20,
                      7 => 25,
                      8 => 25,
                      9 => 50,
                      10 => 50,
                      11 => 50,
                      12 => 60,
                      13 => 60,
                      14 => 60,
                      15 => 70,
                      16 => 100, # higher, tbh
                      17 => 100,
                      18 => 100,
                      19 => 80,
                      20 => 60,
                      21 => 50,
                      22 => 40,
                      23 => 40,
                      24 => 20)

num_nodes = 10000
max_size = 4
hyperedge_colormap = :thermal
linecolors = get(colorschemes[hyperedge_colormap], [1, 2.5, 3], (1, max_size))

# plot every nth point 
skip_middle = Dict(:maj_rts => 1,
                   :maj_rtr => 1,
                   :prop_rts => 5,
                   :prop_rtr => 10)

prompt = true

if plot_type == :total
    filename = "./figures/p_sweep_hyperedges.pdf"
    title = "Stable distribution of hyperedges"
    ylabel = "# of hyperedges"
    ylims = (0, 12400)
    yticks = 0:2500:12000
elseif plot_type == :active
    if plot_analytical
        filename = "./figures/p_sweep_active_analytical.pdf"
    else
        filename = "./figures/p_sweep_active.pdf"
    end
    title = "Parabola height"
    ylabel = L"density of active links $\rho_{d_i}$"
    ylims = (0, 0.5)
    yticks = 0:0.1:1
end

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
xticks = vcat(0.0:0.2:1.0, [0.93])
xlabels = [i != 1 ? "$i" : "$(Int64(i))" for i in xticks]
ax.xticks = (xticks, xlabels)
ax.yticks = yticks
ax.xlabel = L"p"
ax.ylabel = ylabel
xlims!(ax, (-0.03, 1.03))
ylims!(ax, ylims)

p_sweep = vcat(collect(0.0:0.2:0.6), collect(0.65:0.05:0.8), collect(0.85:0.01:1.0))
total_sweep = Dict(size => Float64[] for size in 2:max_size)
total_error = Dict(size => Float64[] for size in 2:max_size)

for (batch_folder, batch_num) in data_folder
    # we first collect the average hyperedge count for each run separately and then take the 
    # average across all runs
    run_values = Dict(size => [] for size in 2:max_size)

    all_magnetization = Dict(size => [] for size in 2:max_size)
    all_data = Dict(size => [] for size in 2:max_size)

    for (run_folder, run_num) in batch_folder
        time, state_count_A = load_data("state_count_A.csv", run_folder)

        for size in 2:max_size
            if plot_type == :total
                indices, values = load_data("hyperedge_count_$size.csv", run_folder)
                # filter to get only the values after the system has converged to a stable state
                values = values[indices .>= starting_point[batch_num]]
                avg = mean(values)
                push!(run_values[size], avg)
            elseif plot_type == :active
                indices, values = load_data("active_hyperedge_count_$size.csv", run_folder)
                values = values ./ num_nodes
                magnetization = 2 .* state_count_A ./ num_nodes .- 1

                # filter to get only the values after the system has converged to a stable state
                has_converged = indices .>= starting_point[batch_num]
                append!(all_magnetization[size], magnetization[has_converged])
                append!(all_data[size], values[has_converged])
            end
        end
    end

    # average over all runs
    for size in 2:max_size
        if plot_type == :total
            average = mean(run_values[size])
            error = std(run_values[size])
            push!(total_sweep[size], average)
            push!(total_error[size], error)
        elseif plot_type == :active
            p, cov = fit_parabola(all_magnetization[size], all_data[size], false)
            peak = p(0)
            error = sqrt(cov[1, 1])
            push!(total_sweep[size], peak)
            push!(total_error[size], error)
        end
    end
end

for (i, size) in enumerate(2:max_size)
    scatter!(ax, p_sweep, total_sweep[size];
             color=linecolors[i], label="size $size", markersize=8, marker=:xcross)
    if draw_error_bars
        errorbars!(ax, p_sweep, total_sweep[size], total_error[size];
                   color=linecolors[i], whiskerwidth=7, linewidth=0.5)
    end
end

# plot the analytical results
if plot_analytical
    params = load_params(joinpath(input_folder, "batch_001", "input_params.json"))
    tspan = (0.0, 500.0)
    p_sweep_analytical = vcat([0, 0.01, 0.02, 0.03], p_sweep[2:end])
    parabola_height_analytical = Dict(s => Float64[] for s in 2:max_size)
    num_total_hyperedges = Dict(s => zero(p_sweep_analytical) for s in 2:max_size)
    for (i, p) in enumerate(p_sweep_analytical)
        params.model_params.adaptivity_prob = p
        t, sol = moment_expansion(params, tspan, moment_closure)

        num_active_hyperedges = Dict(s => 0.0 for s in 2:max_size)
        for label in all_labels(max_size)
            if AH.order(label) != 1
                continue
            end
            size = AH.size(label)[1]
            num_total_hyperedges[size][i] += sol[label][end]
            if label.left[AH.A] == 0 || label.left[AH.B] == 0
                continue
            end
            num_active_hyperedges[size] += sol[label][end]
        end
        for size in 2:max_size
            height = num_active_hyperedges[size] / num_nodes
            push!(parabola_height_analytical[size], height)
        end
    end

    if plot_type == :active
        analytical_solution = parabola_height_analytical
    else
        analytical_solution = num_active_hyperedges
    end
    for size in 2:max_size
        if size == 2
            lines!(ax, p_sweep_analytical, analytical_solution[size];
                   color=linecolors[size - 1],
                   linewidth=1.7,
                   linestyle=:dash, label="mean-field")
        else
            lines!(ax, p_sweep_analytical, analytical_solution[size];
                   color=linecolors[size - 1],
                   linewidth=1.7,
                   linestyle=:dash)
        end
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