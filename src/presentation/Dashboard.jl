export Dashboard, run!, record!, reset!

using DrWatson

#! format: off
"""
A list of dependencies of Panels on Measurements. 
"""
PANEL_DEPENDENCIES = Dict{DataType,
                          Vector{DataType}}(
        HypergraphPanel             => [],
        StateDistPanel              => [StateCount],
        HyperedgeDistPanel          => [HyperedgeCount],
        ActiveHyperedgeDistPanel    => [ActiveHyperedgeCount],
        SlowManifoldPanel           => [StateCount, ActiveHyperedgeCount],
        ActiveLifetimePanel         => [ActiveLifetime],
        FinalMagnetizationPanel     => [FinalMagnetization])
#! format: on

struct Dashboard
    fig::Figure
    panels::Vector{AbstractPanel}
    mo::ModelObservable
    is_interactive::Bool
end

#! format: off
"""
    function Dashboard(model::AbstractModel;
                       panel_types=[StateDistPanel,
                                   HyperedgeDistPanel,
                                   ActiveHyperedgeDistPanel,
                                   SlowManifoldPanel,
                                   ActiveLifetimePanel,
                                   FinalMagnetizationPanel],
                       vparams::VisualizationParams,
                       is_interactive::Bool=false)

A visualization of the evolution of the hypergraph during the simulation.

# Arguments
- `model::AbstractModel` - the model on which the dashboard is based on. 
- `panel_types::Vector{DataType}` - a list of types of `Panel`s to plot. 
- `vparams::VisualizationParams` - parameters used for visualization
- `is_interactive::Bool` - if true, the Dashboard will have interactive GUI elements to controll the simulation. Not implemented yet. 
"""
function Dashboard(model::AbstractModel;
                   panel_types=[StateDistPanel,
                                HyperedgeDistPanel,
                                ActiveHyperedgeDistPanel,
                                SlowManifoldPanel,
                                ActiveLifetimePanel,
                                FinalMagnetizationPanel],
                   vparams::VisualizationParams,
                   is_interactive::Bool=false)
    #! format: on
    fig = Figure(; resolution=(1200, 800))
    display(fig)

    # collect a list of the measurements on which the panels depend on
    measurements = Vector{DataType}()
    for panel in panel_types
        append!(measurements, PANEL_DEPENDENCIES[panel])
    end
    measurement_types = unique(measurements)

    mo = ModelObservable{typeof(model)}(model, measurement_types;
                                        skip_points=vparams.skip_points,
                                        buffer_size=vparams.buffer_size)

    # the plots are in the left half of the figure, the interactive controlls on the right
    plot_box = fig[1, 1] = GridLayout()
    if is_interactive
        controls_box = fig[2, 1] = GridLayout()
    end

    # determine the number of rows and columns 
    num_panels = length(panel_types)
    nrows = Int64(floor(sqrt(num_panels)))
    ncols = Int64(ceil(sqrt(num_panels)))

    panels = []
    graph_properties = Dict(:num_nodes => get_num_nodes(mo.network[]),
                            :max_hyperedge_size => get_max_hyperedge_size(mo.network[]),
                            :num_hyperedges => get_num_hyperedges(mo.network[]))

    for (i, type) in enumerate(panel_types)
        col = mod1(i, ncols)
        row = (i - 1) ÷ ncols + 1

        # HypergraphPanel doesn't take a Measurement and instead needs the network, so 
        # we treat it separately
        if type <: HypergraphPanel
            panel = HypergraphPanel(plot_box[col, row], mo.network,
                                    graph_properties,
                                    vparams)
        else
            # add the dependent measurements to the arg list
            measurement_types = PANEL_DEPENDENCIES[type]
            measurements = Dict()
            for type in measurement_types
                sym = Symbol(_snake_case("$type"))
                measurements[sym] = getfield(mo, sym)
            end
            panel = type(plot_box[col, row], measurements, graph_properties, vparams)
        end
        push!(panels, panel)
    end

    return Dashboard(fig, panels, mo, is_interactive)
end

"""
    run!(dashboard::Dashboard, num_steps::Integer)

Run the simulation for `num_steps` time steps or until the hypergraph runs out of active hyperedges.
"""
function run!(dashboard::Dashboard, num_steps::Integer)
    mo = dashboard.mo

    if dashboard.is_interactive
        # TODO: something should happen here
    else
        for panel in dashboard.panels
            if typeof(panel) <: AbstractTimeSeriesPanel
                panel.xhigh = num_steps
            elseif typeof(panel) <: ActiveLifetimePanel
                panel.yhigh = num_steps^1.05
            end
            set_lims!(panel)
        end

        active_lifetime = num_steps
        for i in 1:num_steps
            step!(mo)
            num_active_hyperedges = get_num_active_hyperedges(mo.network[])
            if i % mo.buffer_size == 0 || num_active_hyperedges == 0
                for panel in dashboard.panels
                    if !(typeof(panel) <: ActiveLifetimePanel)
                        set_lims!(panel)
                    end
                end
            end

            # stop the simulation early if we run out of active hyperdeges
            if num_active_hyperedges == 0
                flush_buffers!(mo)
                active_lifetime = i
                break
            end
        end
        record_measurements!(mo, :run)
    end
    return dashboard
end

"""
    record!(dashboard::Dashboard, filename::String, num_steps::Int64, framerate::Int64)

Run the simulation for `num_steps` time steps and record a video of the dashboard. 

The video is saved to ./videos in a .mp4 format. `filename` should only contain the name 
of the file without the extension.
"""
function record!(dashboard::Dashboard, filename::String, num_steps::Int64,
                 framerate::Int64)
    savepath = joinpath("videos", filename * ".mp4")

    num_updates = num_steps ÷ mo.buffer_size
    record(dashboard.fig, savepath, 1:num_updates; framerate=framerate, compression=1) do i
        return run!(dashboard, mo.buffer_size)
    end
    return dashboard
end

"""
    reset!(dashboard::Dashboard, model::AbstractModel)

Reset the dashboard to run the next simulation from the batch. 

The old history plot lines are made inactive and are grayed out. 
The observables in the ModelObservable are reset to track the new data from `model`.
"""
function reset!(dashboard::Dashboard, model::AbstractModel)

    # gray out the history plot lines
    for panel in dashboard.panels
        if typeof(panel) <: AbstractTimeSeriesPanel ||
           typeof(panel) <: SlowManifoldPanel
            deactivate_lines!(panel)
        end
    end

    # reset observables
    rebind_model!(dashboard.mo, model)
    return dashboard
end

"""
    save(dashboard::Dashboard, folder::String, filename::String)

Save the state of the dashboard as a figure.
The filename should be given with an extension, for example, dash.png. 
"""
function save(dashboard::Dashboard, folder::String, filename::String)
    return GLMakie.save(joinpath(folder, filename), dashboard.fig)
end