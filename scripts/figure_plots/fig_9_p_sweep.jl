using StatsBase

using AdaptiveHypergraphs
const AH = AdaptiveHypergraphs

include("figure_plotting_tools.jl")

# ======================================================
# -------------------- INPUT ---------------------------

draw_error_bars = true

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
              :prop_rts => "Proportional voting\nrewire-to-same",
              :prop_rtr => "Proportional voting\nrewire-to-random")

input_folders = Dict(:maj_rts => folder_maj_rts,
                     :maj_rtr => folder_maj_rtr,
                     :prop_rts => folder_prop_rts,
                     :prop_rtr => folder_prop_rtr)

num_nodes = 10000
max_size = 4
linecolors = hyperedge_linecolors(max_size)

prompt = false

# number of active hyperedges below which the simulation counts as converged
converged_threshold = 10

filename = "./figures/p_sweep.pdf"

xlabel = L"$p$"
ylabel = L"$|m(t \rightarrow \infty)|$"
xticks = 0:0.2:1
yticks = 0:0.2:1
xlims = (-0.02, 1.02)
ylims = (-0.02, 1.2)

max_Δp = 1
max_Δm = 0.06

cache_file = "scripts/figure_plots/cache/analytical_stable_points"

# ------------------------------------------------------
# ======================================================

"""
Compute the absolute magnetization in the absorbing state 
from the simulated data
"""
function simulated_abs_magnetization(batch_folder)
    # first, collect the final magnetization for each run separately, then take the 
    # average across all runs
    magnetization = Float64[]
    params = load_params(joinpath(batch_folder.folder, "input_params.json"))
    p = params.model_params.adaptivity_prob

    for (run_folder, _) in batch_folder
        # check if the simulation has converged (number of acitve hyperedges is close to 0)
        _, values = load_data("active_hyperedge_count_2.csv", run_folder)
        if values[end] > converged_threshold
            continue
        end

        _, state_count_A = load_data("state_count_A.csv", run_folder)
        mag = 2 .* state_count_A ./ num_nodes .- 1
        append!(magnetization, abs(mag[end]))
    end

    return p, mean(magnetization), std(magnetization)
end

"""
Compute the analytical stable point in the (ρ, m)-coordinates 
(active hyperedge density and magnetization). 
"""
function analytical_stable_point(params)
    tspan = (0.0, 200.0)
    num_nodes = params.network_params.num_nodes

    t, sol = moment_expansion(params, tspan, moment_closure)

    num_A_nodes = sol[OrderZeroMotif(AH.A)][end]
    mag = 2 * num_A_nodes / num_nodes - 1

    # get the total number of active hyperedges
    active_hyperedges = 0
    for a in 1:(max_size - 1), b in 1:(max_size - a)
        active_hyperedges += sol[OrderOneMotif(a, b)][end]
    end

    mag = 2 * num_A_nodes / num_nodes - 1
    active_hyperedge_density = active_hyperedges / num_nodes

    return active_hyperedge_density, mag
end

"""
Compute the absolute magnetization in the absorbing state analytically
for majority voting.
    
Here, we just let the system evolve until it reaches an absorbing stable state. 
"""
function analytical_abs_magnetization_maj(params)
    params.network_params.state_A_prob = 0.5001

    cached_values = Dict{Float64,Float64}()
    prop_rule = typeof(params.model_params.adaptivity_rule) <: RewireToRandom ? "rtr" :
                "rts"
    cache_file_local = "$(cache_file)_maj_$(prop_rule).csv"
    if isfile(cache_file_local)
        open(cache_file_local) do io
            for line in eachline(io)
                pair = split(line, ",")
                key = parse(Float64, pair[1])
                value = parse(Float64, pair[2])
                cached_values[key] = value
            end
        end
    end

    # since the magnetization can change very rapidly, we choose 
    # the values of p adaptively. 
    function compute_mag(p)
        if p in keys(cached_values)
            return cached_values[p]
        end
        params.model_params.adaptivity_prob = p
        _, mag = analytical_stable_point(params)
        cached_values[p] = mag
        mkpath(splitdir(cache_file_local)[1])
        open(cache_file_local, "a") do io
            return write(io, "$(p),$(mag)\n")
        end
        return mag
    end

    p, mag = adaptive_sweep((0.0, 1.0), compute_mag, max_Δp, max_Δm)

    p = Vector{Float64}(p)
    mag = Vector{Float64}(mag)

    return p, mag
end

"""
Compute the absolute magnetization in the absorbing state analytically
for proportional voting. 

Here, we calculate multiple points on the slow manifold and fit a parabola 
to the points. The absolute magnetization is estimated from the roots 
of the parabola. 
"""
function analytical_abs_magnetization_prop(params)
    cached_values = Dict{Float64,Float64}()
    adaptivity_rule = params.model_params.adaptivity_rule
    prop_rule = typeof(adaptivity_rule) <: RewireToRandom ? "rtr" : "rts"
    rule_cache_file = "$(cache_file)_prop_$(prop_rule).csv"
    if isfile(rule_cache_file)
        open(rule_cache_file) do io
            for line in eachline(io)
                pair = split(line, ",")
                key = parse(Float64, pair[1])
                value = parse(Float64, pair[2])
                cached_values[key] = value
            end
        end
    end

    # since the magnetization can change very rapidly, we choose 
    # the values of p adaptively. 
    function compute_absorbing_mag(p)
        if p in keys(cached_values)
            return cached_values[p]
        end

        params.model_params.adaptivity_prob = p

        # compute the stable state for different initial magnetizations
        # and fit a parabola to the data points
        magnetization_vec = Float64[]
        active_hyperedges_vec = Float64[]

        # first, compute the height of the parabola at m0 = 0 to see if the parabola exists
        params.network_params.state_A_prob = 0.5
        active_hyperedge_density, mag = analytical_stable_point(params)
        # if the height of the parabola is equal to zero at m0 = 0, then 
        # the magnetization in the absorbing state is also equal to zero. 
        if active_hyperedge_density < 1e-6
            cached_values[p] = 0
            open(rule_cache_file, "a") do io
                return write(io, "$(p),0\n")
            end
            return 0ö
        end

        # if the parabola does exist, compute the magnetization at neighboring points
        # and fit a parabola to the results to obtain the roots. 
        for initial_mag in -0.95:0.05:0.95
            params.network_params.state_A_prob = initial_mag * 0.5 + 0.5

            active_hyperedge_density, mag = analytical_stable_point(params)

            # if the stable state has a non-zero denisty of active hyperedges, use it to compute 
            # the parabola. 
            if active_hyperedge_density > 1e-6
                push!(magnetization_vec, mag)
                push!(active_hyperedges_vec, active_hyperedge_density)
            end
        end

        if length(magnetization_vec) < 3
            mag = 0
            cached_values[p] = mag
            open(rule_cache_file, "a") do io
                return write(io, "$(p),$(mag)\n")
            end
            return mag
        end

        polynomial, cov = fit_parabola(magnetization_vec, active_hyperedges_vec, false)
        @show polynomial

        mag = roots(polynomial)[2]
        cached_values[p] = mag
        mkpath(splitdir(rule_cache_file)[1])
        open(rule_cache_file, "a") do io
            return write(io, "$(p),$(mag)\n")
        end

        return mag
    end

    p, mag = adaptive_sweep((0.0, 1.0), compute_absorbing_mag, max_Δp, max_Δm)

    p = Vector{Float64}(p)
    mag = Vector{Float64}(mag)

    return p, mag
end

# create figure
# fig = create_figure(:large, 2 / 1, 0.8)
fig = create_figure(:large, 3 / 2)
ax_gridpos = fig[1, 1] = GridLayout()

axis = add_four_rules_axes!(ax_gridpos, titles,
                            xlabel, ylabel,
                            xticks, yticks,
                            xlims, ylims)

for panel in panels
    @show panel
    p_sweep = Float64[]
    magnetization = Float64[]
    magnetization_error = Float64[]

    data_folder = DataFolder(input_folders[panel], :simulation)
    for (batch_folder, batch_num) in data_folder
        p, mag, error = simulated_abs_magnetization(batch_folder)
        push!(p_sweep, p)
        push!(magnetization, mag)
        push!(magnetization_error, error)
    end

    scatter!(axis[panel], p_sweep, magnetization;
             markersize=8, marker=:xcross, label="simulation")
    if draw_error_bars
        errorbars!(axis[panel], p_sweep, magnetization, magnetization_error;
                   whiskerwidth=7, linewidth=0.5)
    end

    params = load_params(joinpath(input_folders[panel], "batch_001", "input_params.json"))
    if panel == :maj_rtr || panel == :maj_rts
        p, mag = analytical_abs_magnetization_maj(params)
    else
        p, mag = analytical_abs_magnetization_prop(params)
    end
    lines!(axis[panel], p, mag;
           linewidth=1.7,
           linestyle=:dash, label="mean-field")
end

Legend(fig[1, 2], axis[:prop_rtr]; merge=true, orientation=:vertical, labelsize=12)

display(fig)

if prompt
    prompt_for_save(filename, fig)
end