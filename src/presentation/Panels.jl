export AbstractPanel, AbstractTimeSeriesPanel, HypergraphPanel, StateDistPanel,
       HyperedgeDistPanel, ActiveHyperedgesPanel, SlowManifoldPanel, deactivate_lines!,
       set_lims!

"""
AbstractPanel

A subfigure in the dashboard that visualizes specific data.

In contrast to standard subfigures, each Panel type represents only one specific plot. 
This makes it easier to tweak some plot-specific parameters like axis limits.

Panels often visualize AbstractTimeSeries (see ModelObservable.jl), but there is not a one-to-one 
correspondence: a Panel can visualize multiple AbstractTimeSeries and the same AbstractTimeSeries can 
be used in multiple Panels. 
"""
abstract type AbstractPanel end

"""
AbstractTimeSeriesPanel <: AbstractPanel

A subtype of panel that has time on the x-axis.
"""
abstract type AbstractTimeSeriesPanel <: AbstractPanel end

mutable struct HypergraphPanel <: AbstractPanel
    mo::ModelObservable
    network::HyperNetwork
    axes::Axis
    xlow::Union{Real,Nothing}
    xhigh::Union{Real,Nothing}
    ylow::Union{Real,Nothing}
    yhigh::Union{Real,Nothing}
end

function HypergraphPanel(box::GridSubposition, mo::ModelObservable;
                         node_colormap=:RdYlGn_6, hyperedge_colormap=:thermal,
                         xlow=nothing, xhigh=nothing,
                         ylow=nothing, yhigh=nothing)
    ax, _ = hypergraphplot(box, mo.network; node_colormap,
                           hyperedge_colormap)
    ax.title = "Visualization of the hypergraph"
    return HypergraphPanel(mo, network, ax, xlow, xhigh, ylow, yhigh)
end

mutable struct StateDistPanel <: AbstractTimeSeriesPanel
    time_series::Vector{StateCount}
    axes::Axis
    lines::Vector{Lines}
    xlow::Union{Real,Nothing}
    xhigh::Union{Real,Nothing}
    ylow::Union{Real,Nothing}
    yhigh::Union{Real,Nothing}
end

function StateDistPanel(box::GridSubposition, mo::ModelObservable;
                        node_colormap=:RdYlGn_6,
                        xlow=0, xhigh=nothing,
                        ylow=-10, yhigh=nothing)
    lines = []
    title = "Distribution of states"
    ax = Axis(box; title=title)
    num_states = length(instances(State))
    linecolors = get(colorschemes[node_colormap], 1:num_states, (1, num_states))
    for (i, series) in enumerate(mo.state_series)
        l = lines!(ax,
                   series.observable;
                   label="# of $(series.state) nodes",
                   color=linecolors[i])
        push!(lines, l)
    end
    xlims!(ax; low=xlow, high=xhigh)
    ylims!(ax; low=ylow, high=yhigh)

    axislegend(ax; labelsize=12)

    return StateDistPanel(mo.state_series, ax, lines, xlow, xhigh, ylow, yhigh)
end

mutable struct HyperedgeDistPanel <: AbstractTimeSeriesPanel
    time_series::Vector{HyperedgeCount}
    axes::Axis
    lines::Vector{Lines}
    xlow::Union{Real,Nothing}
    xhigh::Union{Real,Nothing}
    ylow::Union{Real,Nothing}
    yhigh::Union{Real,Nothing}
end

function HyperedgeDistPanel(box::GridSubposition, mo::ModelObservable;
                            hyperedge_colormap=:thermal,
                            xlow=0, xhigh=nothing,
                            ylow=-10, yhigh=nothing)
    lines = []
    title = "Distribution of hyperdeges"
    ax = Axis(box; title=title)
    max_size = get_max_hyperedge_size(mo.network[])
    linecolors = get(colorschemes[hyperedge_colormap], 1:max_size, (1, max_size))
    for (i, series) in enumerate(mo.hyperedge_series)
        l = lines!(ax,
                   series.observable;
                   label="hyperedges of size $(series.size)",
                   color=linecolors[series.size - 1])
        push!(lines, l)
    end
    xlims!(ax; low=xlow, high=xhigh)
    ylims!(ax; low=ylow, high=yhigh)
    axislegend(ax; labelsize=12)

    return HyperedgeDistPanel(mo.hyperedge_series, ax, lines, xlow, xhigh, ylow, yhigh)
end

mutable struct ActiveHyperedgesPanel <: AbstractTimeSeriesPanel
    time_series::Vector{ActiveHyperedgeCount}
    axes::Axis
    lines::Vector{Lines}
    xlow::Union{Real,Nothing}
    xhigh::Union{Real,Nothing}
    ylow::Union{Real,Nothing}
    yhigh::Union{Real,Nothing}
end

function ActiveHyperedgesPanel(box::GridSubposition, mo::ModelObservable;
                               xlow=0, xhigh=nothing,
                               ylow=-10, yhigh=nothing)
    lines = []
    title = "Number of active hyperedges"
    ax = Axis(box[1, 1]; title=title)
    l = lines!(ax,
               mo.active_hyperedges_series.observable)
    push!(lines, l)
    xlims!(ax; low=xlow, high=xhigh)
    ylims!(ax; low=ylow, high=yhigh)
    return ActiveHyperedgesPanel([mo.active_hyperedges_series], ax, lines, xlow, xhigh,
                                 ylow, yhigh)
end

mutable struct SlowManifoldPanel <: AbstractTimeSeriesPanel
    time_series::Vector{AbstractTimeSeries}
    axes::Axis
    lines::Vector{Lines}
    xlow::Union{Real,Nothing}
    xhigh::Union{Real,Nothing}
    ylow::Union{Real,Nothing}
    yhigh::Union{Real,Nothing}
end

function SlowManifoldPanel(box::GridPosition, mo::ModelObservable;
                           xlow=nothing, xhigh=nothing,
                           ylow=nothing, yhigh=nothing)
    lines = []
    title = "Slow manifold plot"
    ax = Axis(box[1, 1]; title=title)
    l = lines!(ax,
               mo.state_series[1].observable,
               mo.active_hyperedges_series.observable)
    xlims!(ax; low=xlow, high=xhigh)
    ylims!(ax; low=ylow, high=yhigh)
    push!(lines, l)
    return SlowManifoldPanel([mo.state_series[1], mo.active_hyperedges_series], ax, lines,
                             xlow, xhigh, ylow, yhigh)
end

mutable struct ActiveLifetimePanel <: AbstractPanel
    run_series::ActiveLifetime
    axes::Axis
    lines::Vector{Scatter}
    xlow::Union{Real,Nothing}
    xhigh::Union{Real,Nothing}
    ylow::Union{Real,Nothing}
    yhigh::Union{Real,Nothing}
end

function ActiveLifetimePanel(box::GridPosition, mo::ModelObservable;
                             xlow=-0.2, xhigh=nothing,
                             ylow=10, yhigh=nothing)
    lines = []
    title = "Time to depletion of active hyperdeges"
    ax = Axis(box[1, 1]; title=title, yscale=log10)
    l = scatter!(ax,
                 mo.active_lifetime.observable)
    xlims!(ax; low=xlow, high=xhigh)
    ylims!(ax; low=ylow, high=yhigh)
    push!(lines, l)
    return ActiveLifetimePanel(mo.active_lifetime, ax, lines, xlow, xhigh, ylow, yhigh)
end

function deactivate_lines!(panel::SlowManifoldPanel)
    lines!(panel.axes,
           panel.time_series[1].observable[],
           panel.time_series[2].observable[];
           linewidth=1,
           color=(:gray, 0.5))
    # Bring the lines tied to observables in front of the gray lines
    for line in panel.lines
        translate!(line, 0, 0, 1)
    end
    return panel
end

function deactivate_lines!(panel::AbstractTimeSeriesPanel)
    for series in panel.time_series
        lines!(panel.axes,
               series.observable[];
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