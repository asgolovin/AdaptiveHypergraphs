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

    xlow, xhigh = (0, nothing)
    ylow, yhigh = (-0.05num_nodes, 1.05num_nodes)

    lines = []
    title = "Distribution of states"
    ax = Axis(box; title=title)
    num_states = length(instances(State))
    linecolors = get(colorschemes[node_colormap], 1:num_states, (1, num_states))

    logs = MeasurementLog[]
    for (i, measurement) in enumerate(state_count)
        log = measurement.log
        l = lines!(ax,
                   log.indices, log.values;
                   label="# of $(measurement.label) nodes",
                   color=linecolors[i])
        push!(lines, l)
        push!(logs, log)
    end
    xlims!(ax; low=xlow, high=xhigh)
    ylims!(ax; low=ylow, high=yhigh)

    axislegend(ax; labelsize=12)

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
    lines = []
    title = "Distribution of hyperdeges"
    ax = Axis(box; title=title)
    linecolors = get(colorschemes[hyperedge_colormap], 1:max_size, (1, max_size))
    logs = MeasurementLog[]
    for (i, measurement) in enumerate(hyperedge_count)
        log = measurement.log
        size = measurement.label
        l = lines!(ax,
                   log.indices, log.values;
                   label="hyperedges of size $(size)",
                   color=linecolors[size - 1])
        push!(lines, l)
        push!(logs, log)
    end
    xlims!(ax; low=xlow, high=xhigh)
    ylims!(ax; low=ylow, high=yhigh)
    axislegend(ax; labelsize=12)

    return HyperedgeDistPanel(logs, ax, lines, xlow, xhigh, ylow, yhigh)
end

# mutable struct FinalHyperedgeDistPanel <: AbstractTimeSeriesPanel
#     measurement_logs::Vector{MeasurementLog}
#     axes::Axis
#     lines::Vector{Lines}
#     xlow::Union{Real,Nothing}
#     xhigh::Union{Real,Nothing}
#     ylow::Union{Real,Nothing}
#     yhigh::Union{Real,Nothing}
# end
# 
# function FinalHyperedgeDistPanel(box::GridPosition;
#                                  final_hyperedge_dist::FinalHyperedgeDist,
#                                  max_size::Int64,
#                                  hyperedge_colormap)
#     (xlow, xhigh) = (0, nothing)
#     (ylow, yhigh) = (-10, nothing)
#     lines = []
#     title = "Final distribution of hyperdeges"
#     ax = Axis(box; title=title)
#     linecolors = get(colorschemes[hyperedge_colormap], 1:max_size, (1, max_size))
#     logs = MeasurementLog[]
#     for i in 1:(max_size - 1)
#         l = lines!(ax,
#                    final_hyperedge_dist.values;)
#         #color=linecolors[i])
#         push!(lines, l)
#         push!()
#     end
#     xlims!(ax; low=xlow, high=xhigh)
#     ylims!(ax; low=ylow, high=yhigh)
# 
#     return FinalHyperedgeDistPanel(final_hyperedge_dist, ax, lines, xlow, xhigh, ylow,
#                                    yhigh)
# end

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

    lines = []
    title = "Number of active hyperedges"
    ax = Axis(box[1, 1]; title=title)
    linecolors = get(colorschemes[hyperedge_colormap], 1:max_size, (1, max_size))
    logs = MeasurementLog[]
    for (i, measurement) in enumerate(active_hyperedge_count)
        log = measurement.log
        size = measurement.label
        l = lines!(ax,
                   log.indices, log.values;
                   label="hyperedges of size $(size)",
                   color=linecolors[size - 1])
        push!(lines, l)
        push!(logs, log)
    end
    xlims!(ax; low=xlow, high=xhigh)
    ylims!(ax; low=ylow, high=yhigh)
    return ActiveHyperedgeDistPanel(logs, ax, lines, xlow, xhigh,
                                    ylow, yhigh)
end

mutable struct SlowManifoldPanel <: AbstractPanel
    measurement_logs::Vector{MeasurementLog}
    axes::Vector{Axis}
    lines::Vector{Lines}
    num_subplots::Int64
    num_cols::Int64
    num_rows::Int64
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
    axes = []
    title = "Slow manifold plot"
    num_subplots = max_size - 1
    num_cols = num_subplots <= 3 ? 1 : 2
    num_rows = Int64(ceil(num_subplots / num_cols))
    linecolors = get(colorschemes[hyperedge_colormap],
                     1:(num_subplots + 1), (1, num_subplots + 1))
    logs = MeasurementLog[]
    push!(logs, state_count[1].log)
    for i in 1:num_subplots
        col = mod1(i, num_cols)
        row = (i - 1) ÷ num_cols + 1
        ax = Axis(box[row, col])
        if i == 1
            ax.title = title
        end
        if i != num_subplots
            hidexdecorations!(ax)
        end
        l = lines!(ax,
                   state_count[1].log.values,
                   active_hyperedge_count[i].log.values;
                   color=linecolors[i])
        xlims!(ax; low=xlow, high=xhigh)
        ylims!(ax; low=ylow, high=yhigh)
        push!(lines, l)
        push!(axes, ax)
        push!(logs, active_hyperedge_count[i].log)
    end

    return SlowManifoldPanel(logs,
                             axes,
                             lines,
                             num_subplots, num_cols, num_rows,
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

function deactivate_lines!(panel::SlowManifoldPanel)
    for i in 1:(panel.num_subplots)
        lines!(panel.axes[i],
               panel.measurement_logs[1].values[],
               panel.measurement_logs[i + 1].values[];
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

function set_lims!(panel::SlowManifoldPanel)
    for i in 1:(panel.num_subplots)
        xlims!(panel.axes[i]; low=panel.xlow, high=panel.xhigh)
        ylims!(panel.axes[i]; low=panel.ylow, high=panel.yhigh)
    end
    return panel
end

function set_lims!(panel::AbstractPanel)
    xlims!(panel.axes; low=panel.xlow, high=panel.xhigh)
    ylims!(panel.axes; low=panel.ylow, high=panel.yhigh)
    return panel
end