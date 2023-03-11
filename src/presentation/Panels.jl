using GLMakie
using ColorSchemes

export AbstractPanel, AbstractTimeSeriesPanel, HypergraphPanel, StateDistPanel,
       HyperedgeDistPanel, MomentClosurePanel, ActiveHyperedgeDistPanel, SlowManifoldPanel,
       deactivate_lines!,
       set_lims!

"""
AbstractPanel

A subfigure in the dashboard that visualizes specific data.

In contrast to standard subfigures, each Panel type represents only one specific plot. 
This makes it easier to tweak some plot-specific parameters like axis limits.

Panels often visualize Measurements (see Measurements.jl), but there is not a one-to-one 
correspondence: a Panel can visualize multiple Measurements and the same Measurement can 
be used in multiple Panels. 
"""
abstract type AbstractPanel end

"""
AbstractTimeSeriesPanel <: AbstractPanel

A subtype of panel that has time on the x-axis.
"""
abstract type AbstractTimeSeriesPanel <: AbstractPanel end

mutable struct HypergraphPanel <: AbstractPanel
    network::Observable{HyperNetwork}
    axes::Axis
    xlow::Union{Real,Nothing}
    xhigh::Union{Real,Nothing}
    ylow::Union{Real,Nothing}
    yhigh::Union{Real,Nothing}
end

# don't do anything for panels which don't need the solution
function set_solution(panel::AbstractPanel, t, sol) end

function HypergraphPanel(box::GridPosition,
                         network::Observable{HyperNetwork},
                         graph_properties::Dict,
                         vparams::VisualizationParams)
    node_colormap = vparams.node_colormap
    hyperedge_colormap = vparams.hyperedge_colormap
    ax, _ = hypergraphplot(box, network; node_colormap,
                           hyperedge_colormap)
    ax.title = "Visualization of the hypergraph"
    xlow = 0
    ylow = 0
    return HypergraphPanel(network, ax, xlow, nothing, ylow, nothing)
end

"""
Helper function to create plots for `AbstractTimeSeriesPanel`s. 
"""
function _plot_time_series(box::GridPosition, measurements::Vector{<:AbstractMeasurement},
                           lims;
                           title::String="",
                           linecolors::Union{Vector,Nothing}=nothing,
                           labels::Vector{String}=String[],
                           plot_legend::Bool=true)
    ax = Axis(box[1, 1]; title=title)
    lines = []
    logs = MeasurementLog[]
    for (i, measurement) in enumerate(measurements)
        kwargs = Dict()
        if labels != []
            kwargs[:label] = labels[i]
        end
        if linecolors !== nothing
            kwargs[:color] = linecolors[i]
        end
        log = measurement.log
        l = lines!(ax,
                   log.indices, log.values;
                   kwargs...)
        push!(lines, l)
        push!(logs, log)
    end
    (xlow, xhigh, ylow, yhigh) = lims
    xlims!(ax; low=xlow, high=xhigh)
    ylims!(ax; low=ylow, high=yhigh)
    if length(labels) > 0 && plot_legend
        axislegend(ax; labelsize=10)
    end
    return ax, lines, logs
end

mutable struct StateDistPanel <: AbstractTimeSeriesPanel
    measurement_logs::Vector{MeasurementLog}
    axes::Axis
    lines::Vector{Lines}
    xlow::Union{Real,Nothing}
    xhigh::Union{Real,Nothing}
    ylow::Union{Real,Nothing}
    yhigh::Union{Real,Nothing}
end

function StateDistPanel(box::GridPosition,
                        measurements::Dict,
                        graph_properties::Dict,
                        vparams::VisualizationParams)
    num_nodes = graph_properties[:num_nodes]
    state_count = measurements[:state_count]
    node_colormap = vparams.node_colormap
    num_states = length(keys(state_count))

    xlow, xhigh = (0, nothing)
    ylow, yhigh = (-0.05num_nodes, 1.05num_nodes)
    lims = (xlow, xhigh, ylow, yhigh)

    title = "Number of nodes in every state"
    linecolors = get(colorschemes[node_colormap], 1:num_states, (1, num_states))
    labels = ["# of $(m.label) nodes" for m in state_count]
    ax, lines, logs = _plot_time_series(box, state_count, lims; title, linecolors, labels)
    return StateDistPanel(logs, ax, lines, xlow, xhigh, ylow, yhigh)
end

mutable struct HyperedgeDistPanel <: AbstractTimeSeriesPanel
    measurement_logs::Vector{MeasurementLog}
    axes::Axis
    lines::Vector{Lines}
    xlow::Union{Real,Nothing}
    xhigh::Union{Real,Nothing}
    ylow::Union{Real,Nothing}
    yhigh::Union{Real,Nothing}
    t_sol::Observable{Vector{Float64}}
    u_sol::Dict{Int64,Observable{Vector{Float64}}}
end

function HyperedgeDistPanel(box::GridPosition,
                            measurements::Dict,
                            graph_properties::Dict,
                            vparams::VisualizationParams)
    hyperedge_count = measurements[:hyperedge_count]
    if :avg_hyperedge_count in keys(measurements)
        plot_average = true
        avg_hyperedge_count = measurements[:avg_hyperedge_count]
    else
        plot_average = false
    end

    max_size = graph_properties[:max_hyperedge_size]
    num_hyperedges = graph_properties[:num_hyperedges]
    hyperedge_colormap = vparams.hyperedge_colormap

    xlow, xhigh = (0, nothing)
    ylow, yhigh = (-0.05num_hyperedges, nothing)
    lims = (xlow, xhigh, ylow, yhigh)

    title = "Number of hyperedges"
    linecolors = get(colorschemes[hyperedge_colormap], 1:(max_size - 1), (1, max_size))
    labels = ["hyperedges of size $(m.label)" for m in hyperedge_count]
    ax, lines, logs = _plot_time_series(box, hyperedge_count, lims; title, linecolors,
                                        labels)

    # plot averages over multiple runs
    if plot_average
        for (i, measurement) in enumerate(avg_hyperedge_count)
            total_mean = lift(measurement.values) do values
                if length(values) > 0
                    return values[end]
                else
                    return 0.0
                end
            end

            hlines!(ax, total_mean; color=:gray)
            xpos = Observable(0.0)
            xpos = lift(hyperedge_count[i].indices) do count
                return max(xpos[], 0.1 * maximum(count; init=-100.0))
            end
            label = @lift "$(round($total_mean, digits=1))"
            text!(xpos, total_mean; text=label, textsize=10, offset=(0, 5))
        end
    end

    # plot the analytical solution
    t_sol = Observable(Float64[])
    u_sol = Dict{Int64,Observable{Vector{Float64}}}()
    for (i, size) in enumerate(2:max_size)
        u_sol[size] = Observable(Float64[])
        lines!(ax, t_sol, u_sol[size]; color=linecolors[i], linewidth=2.0)
    end

    return HyperedgeDistPanel(logs, ax, lines, xlow, xhigh, ylow, yhigh, t_sol, u_sol)
end

function set_solution(panel::HyperedgeDistPanel, t, sol)
    empty!(panel.t_sol[])
    append!(panel.t_sol[], t)
    for size in keys(panel.u_sol)
        empty!(panel.u_sol[size][])
    end
    for motif in keys(sol)
        if order(motif) == 0
            continue
        end
        s = size(motif)
        if length(panel.u_sol[s][]) == 0
            panel.u_sol[s][] = copy(sol[motif])
        else
            panel.u_sol[s][] += sol[motif]
        end
    end
    return panel
end

mutable struct ActiveHyperedgeDistPanel <: AbstractTimeSeriesPanel
    measurement_logs::Vector{MeasurementLog}
    axes::Axis
    lines::Vector{Lines}
    xlow::Union{Real,Nothing}
    xhigh::Union{Real,Nothing}
    ylow::Union{Real,Nothing}
    yhigh::Union{Real,Nothing}
    t_sol::Observable{Vector{Float64}}
    u_sol::Dict{Int64,Observable{Vector{Float64}}}
end

function ActiveHyperedgeDistPanel(box::GridPosition,
                                  measurements::Dict,
                                  graph_properties::Dict,
                                  vparams::VisualizationParams)
    active_hyperedge_count = measurements[:active_hyperedge_count]
    max_size = graph_properties[:max_hyperedge_size]
    num_hyperedges = graph_properties[:num_hyperedges]
    hyperedge_colormap = vparams.hyperedge_colormap

    xlow, xhigh = (0, nothing)
    ylow, yhigh = (-0.05num_hyperedges, nothing)
    lims = (xlow, xhigh, ylow, yhigh)

    title = "Number of active hyperedges"
    linecolors = get(colorschemes[hyperedge_colormap], 1:(max_size - 1), (1, max_size))
    labels = ["hyperedges of size $(m.label)" for m in active_hyperedge_count]
    ax, lines, logs = _plot_time_series(box, active_hyperedge_count, lims; title,
                                        linecolors, labels)

    # plot the analytical solution
    t_sol = Observable(Float64[])
    u_sol = Dict{Int64,Observable{Vector{Float64}}}()
    for (i, size) in enumerate(2:max_size)
        u_sol[size] = Observable(Float64[])
        lines!(ax, t_sol, u_sol[size]; color=linecolors[i], linewidth=2.0)
    end

    return ActiveHyperedgeDistPanel(logs, ax, lines, xlow, xhigh,
                                    ylow, yhigh, t_sol, u_sol)
end

function set_solution(panel::ActiveHyperedgeDistPanel, t, sol)
    empty!(panel.t_sol[])
    append!(panel.t_sol[], t)
    for size in keys(panel.u_sol)
        empty!(panel.u_sol[size][])
    end
    for motif in keys(sol)
        if order(motif) == 0
            continue
        end
        s = size(motif)
        if !is_active(motif)
            continue
        end
        if length(panel.u_sol[s][]) == 0
            panel.u_sol[s][] = copy(sol[motif])
        else
            panel.u_sol[s][] += copy(sol[motif])
        end
    end
    return panel
end

mutable struct FirstOrderMotifCountPanel <: AbstractTimeSeriesPanel
    measurement_logs::Vector{MeasurementLog}
    axes::Axis
    lines::Vector{Lines}
    xlow::Union{Real,Nothing}
    xhigh::Union{Real,Nothing}
    ylow::Union{Real,Nothing}
    yhigh::Union{Real,Nothing}
    t_sol::Observable{Vector{Float64}}
    u_sol::Dict{AbstractMotif,Observable{Vector{Float64}}}
end

function FirstOrderMotifCountPanel(box::GridPosition,
                                   measurements::Dict,
                                   graph_properties::Dict,
                                   vparams::VisualizationParams)
    motif_count = measurements[:motif_count]
    max_size = graph_properties[:max_hyperedge_size]
    num_hyperedges = graph_properties[:num_hyperedges]
    node_colormap = vparams.node_colormap

    xlow, xhigh = (0, nothing)
    ylow, yhigh = (-0.05num_hyperedges, nothing)
    lims = (xlow, xhigh, ylow, yhigh)

    title = "Number of first-order motifs"
    first_order_motif_count = filter(x -> order(x.label) == 1, motif_count)
    first_order_motifs = [m.label for m in first_order_motif_count]
    # labels to use in the legend
    labels = ["$label" for label in first_order_motifs]

    linecolors = []
    colorscheme = colorschemes[node_colormap]
    for motif in first_order_motifs
        numA = motif.A
        numB = motif.B
        ratio = numB / (numA + numB)
        linecolor = get(colorscheme, ratio)
        push!(linecolors, linecolor)
    end

    ax, lines, logs = _plot_time_series(box, first_order_motif_count, lims; title,
                                        linecolors, labels)

    # plot the analytical solution
    t_sol = Observable(Float64[])
    u_sol = Dict{AbstractMotif,Observable{Vector{Float64}}}()
    for (i, motif) in enumerate(first_order_motifs)
        u_sol[motif] = Observable(Float64[])
        lines!(ax, t_sol, u_sol[motif]; color=linecolors[i], linewidth=2.0)
    end

    return FirstOrderMotifCountPanel(logs, ax, lines, xlow, xhigh,
                                     ylow, yhigh, t_sol, u_sol)
end

function set_solution(panel::FirstOrderMotifCountPanel, t, sol)
    empty!(panel.t_sol[])
    append!(panel.t_sol[], t)
    for motif in keys(sol)
        if order(motif) != 1
            continue
        end
        panel.u_sol[motif][] = copy(sol[motif])
    end
    return panel
end

mutable struct SecondOrderMotifCountPanel <: AbstractTimeSeriesPanel
    measurement_logs::Vector{Vector{MeasurementLog}}
    axes::Vector{Axis}
    lines::Vector{Lines}
    xlow::Union{Real,Nothing}
    xhigh::Union{Real,Nothing}
    ylow::Union{Real,Nothing}
    yhigh::Union{Real,Nothing}
end

function SecondOrderMotifCountPanel(box::GridPosition,
                                    measurements::Dict,
                                    graph_properties::Dict,
                                    vparams::VisualizationParams)
    motif_count = measurements[:motif_count]
    max_size = graph_properties[:max_hyperedge_size]
    num_hyperedges = graph_properties[:num_hyperedges]
    node_colormap = vparams.node_colormap

    xlow, xhigh = (0, nothing)
    ylow, yhigh = (-0.05num_hyperedges, nothing)

    axes = Axis[]
    lines = Lines[]
    logs = Vector{Vector{MeasurementLog}}()
    colorscheme = colorschemes[node_colormap]

    for size in 2:max_size
        title = "Number of second-order motifs of size $size"
        size_logs = MeasurementLog[]
        ax = Axis(box[size - 1, 1]; title=title)
        if size != max_size
            hidexdecorations!(ax)
        end
        fixed_size_motif_count = filter(x -> order(x.label) == 2, motif_count)
        fixed_size_motif_count = filter(x -> Base.size(x.label.left_motif) == size,
                                        fixed_size_motif_count)
        for (i, meas) in enumerate(fixed_size_motif_count)
            motif = meas.label
            numA = motif.left_motif.A + motif.right.A
            numB = motif.left_motif.B + motif.right.B
            ratio = numB / (numA + numB)
            linecolor = get(colorscheme, ratio)
            l = lines!(ax, meas.indices, meas.values; label="$motif", color=linecolor)
            push!(lines, l)
            push!(size_logs, meas.log)
        end
        push!(axes, ax)
        push!(logs, size_logs)
    end

    return SecondOrderMotifCountPanel(logs, axes, lines, xlow, xhigh,
                                      ylow, yhigh)
end

mutable struct MomentClosurePanel <: AbstractTimeSeriesPanel
    measurement_logs::Vector{Vector{MeasurementLog}}
    axes::Axis
    lines::Vector{Lines}
    xlow::Union{Real,Nothing}
    xhigh::Union{Real,Nothing}
    ylow::Union{Real,Nothing}
    yhigh::Union{Real,Nothing}
end

function MomentClosurePanel(box::GridPosition,
                            measurements::Dict,
                            graph_properties::Dict,
                            vparams::VisualizationParams)
    motif_count = measurements[:motif_count]
    state_count = measurements[:state_count]
    colormap = vparams.misc_colormap

    xlow, xhigh = (0, nothing)
    ylow, yhigh = (0, nothing)

    title = "The true moment closure and the approximation"
    ax = Axis(box; title=title)

    # to reduce the number of options, consider only triples which:
    # - are of order 2, otherwise they are not triples
    # - intersect in an A-node (B-nodes should behave symmetrically)
    triples = filter(x -> order(x.label) == 2, motif_count)
    int_state = A
    filter!(x -> x.label.int[Symbol(int_state)] == 1, triples)

    num_int_nodes = filter(x -> x.label == int_state, state_count)[1].values
    num_triples = length(triples)
    linecolors = get(colorschemes[colormap], 1:num_triples, (1, num_triples))

    lines = Lines[]
    logs = Vector{Vector{MeasurementLog}}()

    for (i, triple) in enumerate(triples)
        motif = triple.label

        # The moment closure is done by replacing [X|Y|Z] by [XY] * [YZ] / [Y]. 
        # Here, we determine XY and YZ and the corresponding labels.
        left_motif = motif.left_motif
        right_motif = motif.right_motif

        # [XY]
        left_motif_values = filter(x -> x.label == left_motif, motif_count)[1].values
        # [YZ]
        right_motif_values = filter(x -> x.label == right_motif, motif_count)[1].values

        # compute the combinatorical prefactor
        int_state = motif.int.A > 0 ? A : B
        if int_state == A
            left_count = left_motif.A
            right_count = right_motif.A
        else
            left_count = left_motif.B
            right_count = right_motif.B
        end
        issymmetrical = left_motif == right_motif ? 0.5 : 1.0
        prefactor = issymmetrical * left_count * right_count

        # The resulting closure approximation
        prediction = @lift prefactor * $left_motif_values .* $right_motif_values ./
                           $num_int_nodes

        l1 = lines!(ax, triple.indices, prediction;
                    label="$left_motif * $right_motif / [ $int_state ]",
                    color=linecolors[i])

        # The actual number of triples
        l2 = lines!(ax, triple.indices, triple.values;
                    label="$motif", color=linecolors[i] * 0.5, linewidth=0.7)

        push!(lines, l1)
        push!(lines, l2)
        triple_log = MeasurementLog[]

        prediction_log = MeasurementLog{Float64,Float64}(triple.indices, prediction)
        push!(triple_log, prediction_log)
        push!(triple_log, triple.log)
        push!(logs, triple_log)
    end

    axislegend()

    return MomentClosurePanel(logs, ax, lines, xlow, xhigh,
                              ylow, yhigh)
end

mutable struct FakeDiffEqPanel <: AbstractTimeSeriesPanel
    measurement_logs::Vector{MeasurementLog}
    axes::Axis
    lines::Vector{Lines}
    xlow::Union{Real,Nothing}
    xhigh::Union{Real,Nothing}
    ylow::Union{Real,Nothing}
    yhigh::Union{Real,Nothing}
end

function FakeDiffEqPanel(box::GridPosition,
                         measurements::Dict,
                         graph_properties::Dict,
                         vparams::VisualizationParams)
    fake_diff_eq = measurements[:fake_diff_eq]
    motif_count = measurements[:motif_count]
    num_hyperedges = graph_properties[:num_hyperedges]
    colormap = vparams.misc_colormap

    plot_conditions = x -> x.label.A > 0
    #plot_conditions = x -> true

    fake_diff_eq = filter(plot_conditions, fake_diff_eq)

    xlow, xhigh = (0, nothing)
    ylow, yhigh = (-0.05num_hyperedges, nothing)
    lims = (xlow, xhigh, ylow, yhigh)

    labels = ["$(m.label)" for m in fake_diff_eq]
    title = "Fake differential equation"

    order_one_motifs = filter(x -> order(x.label) == 1, motif_count)
    filter!(plot_conditions, order_one_motifs)
    num_motifs = length(order_one_motifs)

    linecolors = get(colorschemes[colormap], 1:num_motifs, (1, num_motifs))

    ax, lines, logs = _plot_time_series(box, fake_diff_eq, lims; title,
                                        linecolors, labels, plot_legend=false)

    # plot the true values
    for (i, motif) in enumerate(order_one_motifs)
        l = lines!(ax, motif.indices, motif.values;
                   color=linecolors[i] * 0.6,
                   label="True $(motif.label)",
                   linewidth=0.7)
        push!(lines, l)
        push!(logs, motif.log)
    end

    axislegend()

    return FakeDiffEqPanel(logs, ax, lines, xlow, xhigh,
                           ylow, yhigh)
end

mutable struct ActiveRatioPanel <: AbstractTimeSeriesPanel
    measurement_logs::Vector{MeasurementLog}
    axes::Axis
    lines::Vector{Lines}
    xlow::Union{Real,Nothing}
    xhigh::Union{Real,Nothing}
    ylow::Union{Real,Nothing}
    yhigh::Union{Real,Nothing}
end

function ActiveRatioPanel(box::GridPosition,
                          measurements::Dict,
                          graph_properties::Dict,
                          vparams::VisualizationParams)
    active_hyperedge_count = measurements[:active_hyperedge_count]
    max_size = graph_properties[:max_hyperedge_size]
    num_hyperedges = graph_properties[:num_hyperedges]
    hyperedge_colormap = vparams.hyperedge_colormap

    xlow, xhigh = (0, nothing)
    ylow, yhigh = (-1, nothing)

    title = "Ratio of active hyperedges"
    linecolors = get(colorschemes[hyperedge_colormap], 1:(max_size - 1), (1, max_size))
    labels = ["hyperedges of size $(m.label)" for m in active_hyperedge_count]

    size_two_hyperedges = active_hyperedge_count[1].values

    lines = []
    logs = []

    ax = Axis(box[1, 1]; title=title)
    for (i, measurement) in enumerate(active_hyperedge_count)
        log = measurement.log
        size = measurement.label

        if size == 2
            # for hyperedges of size 2 the ratio is fixed to one
            ratio = lift(x -> ones(length(x)), log.indices)
        else
            ratio = lift((x, y) -> x ./ y, log.values, size_two_hyperedges)
        end
        ratio_log = MeasurementLog{Float64,Float64}(log.indices, ratio)

        l = lines!(ax,
                   log.indices, ratio;
                   label="hyperedges of size $(size)",
                   color=linecolors[size - 1])
        push!(lines, l)
        push!(logs, ratio_log)
    end

    return ActiveRatioPanel(logs, ax, lines, xlow, xhigh, ylow, yhigh)
end

mutable struct SlowManifoldPanel <: AbstractPanel
    measurement_logs::Vector{MeasurementLog}
    axes::Axis
    lines::Vector{Lines}
    xlow::Union{Real,Nothing}
    xhigh::Union{Real,Nothing}
    ylow::Union{Real,Nothing}
    yhigh::Union{Real,Nothing}
end

function SlowManifoldPanel(box::GridPosition,
                           measurements::Dict,
                           graph_properties::Dict,
                           vparams::VisualizationParams)
    state_count = measurements[:state_count]
    active_hyperedge_count = measurements[:active_hyperedge_count]

    if :slow_manifold_fit in keys(measurements)
        plot_fit = true
        slow_manifold_fit = measurements[:slow_manifold_fit]
    else
        plot_fit = false
    end

    max_size = graph_properties[:max_hyperedge_size]
    num_nodes = graph_properties[:num_nodes]
    num_hyperedges = graph_properties[:num_hyperedges]
    hyperedge_colormap = vparams.hyperedge_colormap

    xlow, xhigh = (-0.05num_nodes, 1.05num_nodes)
    ylow, yhigh = (-0.05num_hyperedges, nothing)

    lines = []
    title = "Slow manifold plot"
    linecolors = get(colorschemes[hyperedge_colormap], 1:max_size, (1, max_size))
    logs = MeasurementLog[]
    push!(logs, state_count[1].log)
    ax = Axis(box[1, 1]; title=title)
    for (i, measurement) in enumerate(active_hyperedge_count)
        active_hyperedge_log = measurement.log
        size = measurement.label
        l = lines!(ax,
                   logs[1].values, active_hyperedge_log.values;
                   label="hyperedges of size $(size)",
                   color=linecolors[size - 1])
        push!(lines, l)
        push!(logs, active_hyperedge_log)
    end

    if plot_fit
        for (i, measurement) in enumerate(slow_manifold_fit)
            x = 0:1:num_nodes
            coeffs = lift(measurement.values) do values
                if length(values) == 0
                    return (0.0, 0.0, 0.0)
                else
                    return values[end]
                end
            end
            a = @lift ($coeffs)[1]
            b = @lift ($coeffs)[2]
            c = @lift ($coeffs)[3]
            y = @lift @. $a + $b * x + $c * x^2
            l = lines!(ax, x, y; color=:red)
            translate!(l, 0, 0, 2)
        end
    end

    return SlowManifoldPanel(logs,
                             ax,
                             lines,
                             xlow, xhigh, ylow, yhigh)
end

mutable struct ActiveLifetimePanel <: AbstractPanel
    measurement_logs::Vector{MeasurementLog}
    axes::Axis
    lines::Vector{Scatter}
    xlow::Union{Real,Nothing}
    xhigh::Union{Real,Nothing}
    ylow::Union{Real,Nothing}
    yhigh::Union{Real,Nothing}
end

function ActiveLifetimePanel(box::GridPosition,
                             measurements::Dict,
                             graph_properties::Dict,
                             vparams::VisualizationParams)
    active_lifetime = measurements[:active_lifetime][1]

    xlow, xhigh = (-0.2, nothing)
    ylow, yhigh = (10, nothing)
    lines = []
    title = "Time to depletion of active hyperdeges"
    ax = Axis(box[1, 1]; title=title, yscale=log10)

    l = scatter!(ax,
                 active_lifetime.values)
    xlims!(ax; low=xlow, high=xhigh)
    ylims!(ax; low=ylow, high=yhigh)
    push!(lines, l)
    return ActiveLifetimePanel([active_lifetime.log], ax, lines, xlow, xhigh, ylow, yhigh)
end

mutable struct FinalMagnetizationPanel <: AbstractPanel
    measurement_logs::Vector{MeasurementLog}
    axes::Axis
    lines::Vector{Scatter}
    xlow::Union{Real,Nothing}
    xhigh::Union{Real,Nothing}
    ylow::Union{Real,Nothing}
    yhigh::Union{Real,Nothing}
end

function FinalMagnetizationPanel(box::GridPosition,
                                 measurements::Dict,
                                 graph_properties::Dict,
                                 vparams::VisualizationParams)
    final_magnetization = measurements[:final_magnetization][1]
    num_nodes = graph_properties[:num_nodes]

    xlow, xhigh = (0, nothing)
    ylow, yhigh = (-1.05num_nodes, 1.05num_nodes)
    lines = []
    title = "Final magnetization after a simulation"
    ax = Axis(box[1, 1]; title=title)
    l = scatter!(ax,
                 final_magnetization.values)
    xlims!(ax; low=xlow, high=xhigh)
    ylims!(ax; low=ylow, high=yhigh)
    push!(lines, l)
    return FinalMagnetizationPanel([final_magnetization.log], ax, lines, xlow, xhigh, ylow,
                                   yhigh)
end

mutable struct AvgHyperedgeCountPanel <: AbstractPanel
    measurement_logs::Vector{MeasurementLog}
    axes::Axis
    lines::Vector{Scatter}
    xlow::Union{Real,Nothing}
    xhigh::Union{Real,Nothing}
    ylow::Union{Real,Nothing}
    yhigh::Union{Real,Nothing}
end

function AvgHyperedgeCountPanel(box::GridPosition,
                                measurements::Dict,
                                graph_properties::Dict,
                                vparams::VisualizationParams)
    avg_hyperedge_count = measurements[:avg_hyperedge_count]
    slow_manifold_fit = measurements[:slow_manifold_fit]
    max_size = graph_properties[:max_hyperedge_size]
    hyperedge_colormap = vparams.hyperedge_colormap

    xlow, xhigh = (0, nothing)
    ylow, yhigh = (-0.05max_size, nothing)
    lines = []
    logs = []
    title = "Average number of hyperedges"
    ax = Axis(box[1, 1]; title=title)
    linecolors = get(colorschemes[hyperedge_colormap], 1:(max_size - 1), (1, max_size))

    for size in 2:max_size
        total_count = avg_hyperedge_count[size - 1].values
        active_count = lift(slow_manifold_fit[size - 1].values) do coeffs
            count = Float64[]
            for (a, b, c) in coeffs
                x = -b / (2 * c)
                push!(count, a + b * x + c * x^2)
            end
            return count
        end

        l_total = scatter!(ax,
                           total_count;
                           marker=:circle,
                           color=linecolors[size - 1],
                           markersize=16)
        lines!(ax, total_count; color=linecolors[size - 1])
        l_active = scatter!(ax,
                            active_count;
                            marker=:cross,
                            color=linecolors[size - 1],
                            markersize=16)
        lines!(ax, active_count; color=linecolors[size - 1])
        push!(logs, avg_hyperedge_count[size - 1].log)
        push!(lines, l_total)
        push!(lines, l_active)
    end
    xlims!(ax; low=xlow, high=xhigh)
    ylims!(ax; low=ylow, high=yhigh)
    return AvgHyperedgeCountPanel(logs, ax, lines, xlow, xhigh, ylow,
                                  yhigh)
end

function _bring_to_front!(lines::Vector{Lines})
    # Bring the lines tied to observables in front of the gray lines
    for line in lines
        translate!(line, 0, 0, 1)
    end
end

function deactivate_lines!(panel::SlowManifoldPanel)
    linecolors = [line.attributes.color for line in panel.lines]
    for i in 2:length(panel.measurement_logs)
        lines!(panel.axes,
               panel.measurement_logs[1].values[],
               panel.measurement_logs[i].values[];
               linewidth=1,
               color=(linecolors[i - 1], 0.5))
    end
    _bring_to_front!(panel.lines)
    return panel
end

function deactivate_lines!(panel::SecondOrderMotifCountPanel)
    linecolors = [line.attributes.color for line in panel.lines]
    line_counter = 1
    for (i, ax) in enumerate(panel.axes)
        for log in panel.measurement_logs[i]
            lines!(ax,
                   log.indices[],
                   log.values[];
                   linewidth=1,
                   color=(linecolors[line_counter], 0.5))
            line_counter += 1
        end
    end
    _bring_to_front!(panel.lines)
    return panel
end

function deactivate_lines!(panel::AbstractTimeSeriesPanel)
    if typeof(panel.measurement_logs) <: Vector{Vector{MeasurementLog}}
        logs = reduce(vcat, panel.measurement_logs)
    else
        logs = panel.measurement_logs
    end
    linecolors = [line.attributes.color for line in panel.lines]
    for (i, log) in enumerate(logs)
        lines!(panel.axes,
               log.indices[],
               log.values[];
               linewidth=0.5,
               color=(linecolors[i], 0.5))
    end
    _bring_to_front!(panel.lines)
    return panel
end

function set_lims!(panel::AbstractPanel)
    xlims!(panel.axes; low=panel.xlow, high=panel.xhigh)
    ylims!(panel.axes; low=panel.ylow, high=panel.yhigh)
    return panel
end

function set_lims!(panel::SecondOrderMotifCountPanel)
    for ax in panel.axes
        xlims!(ax; low=panel.xlow, high=panel.xhigh)
        ylims!(ax; low=panel.ylow, high=panel.yhigh)
    end
    return panel
end