export Dashboard, run!, record!, reset!

struct Dashboard
    fig::Figure
    panels::Vector{AbstractPanel}
    mo::ModelObservable
    is_interactive::Bool
end

"""
    Dashboard(model::AbstractModel;
              plot_hypergraph::Bool=false,
              plot_states::Bool=true,
              plot_hyperedges::Bool=true,
              plot_active_hyperedges::Bool=true,
              is_interactive::Bool=false,
              node_colormap = :RdYlGn_6,
              hyperedge_colormap = :thermal)

A visualization of the evolution of the hypergraph during the simulation.
"""
function Dashboard(model::AbstractModel;
                   plot_hypergraph::Bool=false,
                   plot_states::Bool=true,
                   plot_hyperedges::Bool=true,
                   plot_active_hyperedges::Bool=true,
                   plot_slow_manifold::Bool=true,
                   plot_active_lifetime::Bool=true,
                   plot_final_magnetization::Bool=true,
                   is_interactive::Bool=false,
                   node_colormap=:RdYlGn_6,
                   hyperedge_colormap=:thermal)
    fig = Figure(; resolution=(1200, 800))
    display(fig)
    panels = []

    mo = ModelObservable{typeof(model)}(model)

    plot_box = fig[1, 1] = GridLayout()
    if is_interactive
        controls_box = fig[2, 1] = GridLayout()
    end

    # create columns for different panel types
    col_count = 0
    if plot_hypergraph
        col_count += 1
        hg_box = plot_box[1, col_count]
    end
    if plot_states || plot_hyperedges || plot_active_hyperedges
        col_count += 1
        history_box = plot_box[1, col_count]
    end
    if plot_slow_manifold
        col_count += 1
        slow_manifold_box = plot_box[1, col_count]
    end
    if plot_active_lifetime || plot_final_magnetization
        col_count += 1
        run_box = plot_box[1, col_count]
    end

    if plot_hypergraph
        panel = HypergraphPanel(hg_box[1, 1], network; node_colormap, hyperedge_colormap)
        push!(panels, panel)
    end

    if plot_states
        panel = StateDistPanel(history_box[1, 1], mo; node_colormap,
                               ylow=-0.05get_num_nodes(mo.network[]),
                               yhigh=1.05get_num_nodes(mo.network[]))
        push!(panels, panel)
    end

    if plot_hyperedges
        panel = HyperedgeDistPanel(history_box[2, 1], mo; hyperedge_colormap,
                                   ylow=-0.05get_num_hyperedges(mo.network[]))
        push!(panels, panel)
    end

    if plot_active_hyperedges
        active_panel = ActiveHyperedgesPanel(history_box[3, 1], mo; hyperedge_colormap,
                                             ylow=-0.05get_num_active_hyperedges(mo.network[]))
        push!(panels, active_panel)
    end

    if plot_slow_manifold
        max_mag = get_num_nodes(mo.network[])
        panel = SlowManifoldPanel(slow_manifold_box[1, 1], mo;
                                  xlow=-0.05max_mag,
                                  xhigh=1.05max_mag,
                                  ylow=-0.01get_num_hyperedges(mo.network[]))
        push!(panels, panel)
    end

    if plot_active_lifetime
        panel = ActiveLifetimePanel(run_box[1, 1], mo)
        push!(panels, panel)
    end

    if plot_final_magnetization
        panel = FinalMagnetizationPanel(run_box[2, 1], mo;
                                        ylow=-1.05get_num_nodes(mo.network[]),
                                        yhigh=1.05get_num_nodes(mo.network[]))
        push!(panels, panel)
    end

    return Dashboard(fig, panels, mo, is_interactive)
end

"""
    run!(dashboard::Dashboard, num_steps::Integer, steps_per_update::Integer)

Run the simulation for `num_steps` time steps or until the hypergraph runs out of active hyperedges. 
    
The visualization is updated only once every number of steps given by `steps_per_update`.
"""
function run!(dashboard::Dashboard, num_steps::Integer, steps_per_update::Integer)
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
            if i % steps_per_update == 0 || num_active_hyperedges == 0
                flush_buffers!(mo)
                sleep(0.01)
                notify(mo.network)
                for panel in dashboard.panels
                    if !(typeof(panel) <: ActiveLifetimePanel)
                        set_lims!(panel)
                    end
                end
            end

            # stop the simulation early if we run out of active hyperdeges
            if num_active_hyperedges == 0
                active_lifetime = i
                break
            end
        end
        record_active_lifetime!(mo, active_lifetime)
        record_final_magnetization!(mo)
    end
    return dashboard
end

"""
    record!(dashboard::Dashboard, filename::String, num_steps::Integer, steps_per_update::Integer, framerate::Integer)

Run the simulation for `num_steps` time steps and record a video of the dashboard. 

The video is saved to ./videos in a .mp4 format. `filename` should only contain the name 
of the file without the extension.
"""
function record!(dashboard::Dashboard, filename::String, num_steps::Integer,
                 steps_per_update::Integer, framerate::Integer)
    savepath = joinpath("videos", filename * ".mp4")

    num_updates = num_steps ÷ steps_per_update
    record(dashboard.fig, savepath, 1:num_updates; framerate=framerate, compression=1) do i
        return run!(dashboard, steps_per_update, steps_per_update)
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