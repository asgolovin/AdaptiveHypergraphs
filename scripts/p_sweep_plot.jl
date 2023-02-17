using ColorSchemes
using StatsBase

include("./plotting_tools.jl")

input_folder = joinpath(projectdir(),
                        "results/run_2023-02-08_03-12-17_p_sweep")

hyperedge_count = create_measurements(joinpath(input_folder, "batch_001", "run_001"),
                                      HyperedgeCount)
active_hyperedge_count = create_measurements(joinpath(input_folder, "batch_001", "run_001"),
                                             ActiveHyperedgeCount)
state_count = create_measurements(joinpath(input_folder, "batch_001", "run_001"),
                                  StateCount)

params = load_params(joinpath(input_folder, "batch_001", "input_params.json"))

# time from which the trajectories converged to a stable distribution for each batch
# (scientifically measured by looking at the plots)
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

# All points which deviate from the middle of the parabola by at most 0.5 * window_width * num_nodes
# are used to calculate the height of the parabola. 
window_width = 0.10

# Information about the network
max_size = length(params.network_params.num_hyperedges) + 1
num_nodes = params.network_params.num_nodes

# Plotting stuff
fig = Figure()
ax_total_count = Axis(fig[1, 1])
ax_active_count = Axis(fig[2, 1])

p_sweep = Float64[]
total_he_sweep = Dict(size => Float64[] for size in 2:max_size)
total_he_error = Dict(size => Float64[] for size in 2:max_size)
active_he_sweep = Dict(size => Float64[] for size in 2:max_size)
active_he_error = Dict(size => Float64[] for size in 2:max_size)

hyperedge_colormap = params.visualization_params.hyperedge_colormap
linecolors = get(colorschemes[hyperedge_colormap], 1:(max_size - 1), (1, max_size))

for batchdir in readdir(input_folder; join=true)
    if !isdir(batchdir)
        continue
    end
    batch_num = parse(Int64, match(r"batch_([0-9]+)", batchdir)[1])

    params = load_params(joinpath(batchdir, "input_params.json"))

    p = params.model_params.adaptivity_prob
    push!(p_sweep, p)

    # we first collect the average hyperedge count for each run separately and then take the 
    # average across all runs
    run_hyperedge_values = Dict(size => [] for size in 2:max_size)
    run_active_values = Dict(size => [] for size in 2:max_size)

    for rundir in readdir(batchdir; join=true)
        if !isdir(rundir)
            continue
        end
        run_num = parse(Int64, match(r"run_([0-9]+)", splitdir(rundir)[end])[1])

        load_meas(rundir, hyperedge_count)
        load_meas(rundir, active_hyperedge_count)
        load_meas(rundir, state_count)

        for meas in hyperedge_count
            size = meas.label
            indices = meas.indices[]
            values = meas.values[]
            # filter to get only the values after the system has converged to a stable state
            values = values[indices .>= starting_point[batch_num]]
            avg = mean(values)
            push!(run_hyperedge_values[size], avg)
        end

        for meas in active_hyperedge_count
            size = meas.label
            indices = meas.indices[]
            values = meas.values[]
            # mask the values where the magnetization is too far away from the middle
            abs_magnetization = abs.(state_count[1].values[] .- 0.5 * num_nodes) ./
                                num_nodes
            magnetization_mask = abs_magnetization .<= window_width * 0.5

            # filter to get only the values after the system has converged to a stable state
            has_converged = indices .>= starting_point[batch_num]
            values = values[magnetization_mask .& has_converged]

            avg = mean(values)
            if !isnan(avg)
                push!(run_active_values[size], avg)
            end
        end
    end

    # average over all runs
    for size in 2:max_size
        hyperedge_avg = mean(run_hyperedge_values[size])
        hyperedge_std = std(run_hyperedge_values[size])
        active_avg = mean(run_active_values[size])
        active_std = std(run_active_values[size])
        push!(total_he_sweep[size], hyperedge_avg)
        push!(total_he_error[size], hyperedge_std)
        push!(active_he_sweep[size], active_avg)
        push!(active_he_error[size], active_std)
    end
end

for (i, size) in enumerate(2:max_size)
    scatter!(ax_total_count, p_sweep, total_he_sweep[size];
             color=linecolors[i], label="hyperedges of size $size")
    errorbars!(ax_total_count, p_sweep, total_he_sweep[size], total_he_error[size];
               color=linecolors[i], whiskerwidth=10)

    scatter!(ax_active_count, p_sweep, active_he_sweep[size];
             color=linecolors[i], label="hyperedges of size $size")
    errorbars!(ax_active_count, p_sweep, active_he_sweep[size], active_he_error[size];
               color=linecolors[i], whiskerwidth=10)
end

ax_total_count.xlabel = L"p"
ax_total_count.ylabel = "number of hyperedges"
axislegend(ax_total_count)
display(fig)