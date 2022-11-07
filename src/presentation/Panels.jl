export AbstractPanel, AbstractTimeSeriesPanel, HypergraphPanel, StateDistPanel,
       HyperedgeDistPanel, ActiveHyperedgeDistPanel, SlowManifoldPanel, deactivate_lines!,
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
    network::HyperNetwork
    axes::Axis
    xlow::Union{Real,Nothing}
    xhigh::Union{Real,Nothing}
    ylow::Union{Real,Nothing}
    yhigh::Union{Real,Nothing}
end

function HypergraphPanel(box::GridPosition;
                         network::Observable{HyperNetwork},
                         node_colormap=:RdYlGn_6,
                         hyperedge_colormap=:thermal)
    ax, _ = hypergraphplot(box, mo.network; node_colormap,
                           hyperedge_colormap)
    ax.title = "Visualization of the hypergraph"
    return HypergraphPanel(network, ax, xlow, xhigh, ylow, yhigh)
end

"""
Helper function to create plots for `AbstractTimeSeriesPanel`s. 
"""
function _plot_time_series(box::GridPosition, measurements::Vector{<:AbstractMeasurement},
                           lims;
                           title::String="",
                           linecolors::Union{Vector,Nothing}=nothing,
                           labels::Vector{String}=String[])
    ax = Axis(box[1, 1]; title=title)
    lines = []
    logs = MeasurementLog[]
    for (i, measurement) in enumerate(measurements)
        kwargs = Dict()
        if labels != []
            kwargs[:label] = labels[i]
        end
        if linecolors != nothing
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
    if length(labels) > 0
        axislegend(ax; labelsize=12)
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
    ax, lines, logs = _plot_time_series(box, state_count, lims; title, linecolors)
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
end

function HyperedgeDistPanel(box::GridPosition,
                            measurements::Dict,
                            graph_properties::Dict,
                            vparams::VisualizationParams)
    hyperedge_count = measurements[:hyperedge_count]
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

    return HyperedgeDistPanel(logs, ax, lines, xlow, xhigh, ylow, yhigh)
end

mutable struct ActiveHyperedgeDistPanel <: AbstractTimeSeriesPanel
    measurement_logs::Vector{MeasurementLog}
    axes::Axis
    lines::Vector{Lines}
    xlow::Union{Real,Nothing}
    xhigh::Union{Real,Nothing}
    ylow::Union{Real,Nothing}
    yhigh::Union{Real,Nothing}
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

    return ActiveHyperedgeDistPanel(logs, ax, lines, xlow, xhigh,
                                    ylow, yhigh)
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
                 active_lifetime.log.values)
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
                 final_magnetization.log.values)
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
        total_count = @lift [v.total for v in $(avg_hyperedge_count[size - 1].log.values)]
        active_count = @lift [v.active for v in $(avg_hyperedge_count[size - 1].log.values)]
        l_total = scatter!(ax,
                           total_count;
                           marker=:circle,
                           color=linecolors[size - 1])
        l_active = scatter!(ax,
                            active_count;
                            marker=:cross,
                            color=linecolors[size - 1])
        push!(logs, avg_hyperedge_count[size - 1].log)
        push!(lines, l_total)
        push!(lines, l_active)
    end
    xlims!(ax; low=xlow, high=xhigh)
    ylims!(ax; low=ylow, high=yhigh)
    return AvgHyperedgeCountPanel(logs, ax, lines, xlow, xhigh, ylow,
                                  yhigh)
end

function deactivate_lines!(panel::SlowManifoldPanel)
    for i in 2:length(panel.measurement_logs)
        lines!(panel.axes,
               panel.measurement_logs[1].values[],
               panel.measurement_logs[i].values[];
               linewidth=1,
               color=(:gray, 0.5))
    end
    # Bring the lines tied to observables in front of the gray lines
    for line in panel.lines
        translate!(line, 0, 0, 1)
    end
    return panel
end

function deactivate_lines!(panel::AbstractTimeSeriesPanel)
    for log in panel.measurement_logs
        lines!(panel.axes,
               log.indices[],
               log.values[];
               linewidth=1,
               color=(:gray, 0.5))
    end
    # Bring the lines tied to observables in front of the gray lines
    for line in panel.lines
        translate!(line, 0, 0, 1)
    end
    return panel
end

function set_lims!(panel::AbstractPanel)
    xlims!(panel.axes; low=panel.xlow, high=panel.xhigh)
    ylims!(panel.axes; low=panel.ylow, high=panel.yhigh)
    return panel
end