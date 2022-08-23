export Panel, Dashboard, run!, record!, reset!

"""
Types of data that can be visualized in the dashboard.
"""
@enum Panel begin
    hypergraphPanel
    stateDistPanel
    hyperedgeDistPanel
    activeHyperedgesPanel
end


struct Dashboard
    fig::Figure
    panels::Vector{Panel}
    axes::Dict{Panel, Axis}
    lines::Dict{Panel, Any}
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
                   is_interactive::Bool=false,
                   node_colormap = :RdYlGn_6,
                   hyperedge_colormap = :thermal)

    fig = Figure(resolution = (1000, 600))
    display(fig)
    axes = Dict{Panel, Axis}()
    panels = Panel[]
    obs_lines = Dict{Panel, Any}()

    mo = ModelObservable{typeof(model)}(model)

    plot_box = fig[1, 1] = GridLayout()
    if is_interactive
        controls_box = fig[2, 1] = GridLayout()
    end

    plot_count = 0
    
    if plot_hypergraph
        plot_count += 1
        hg_box = plot_box[1, plot_count]
        hgax, _ = hypergraphplot(hg_box[1, 1], mo.network; node_colormap, hyperedge_colormap)
        push!(panels, hypergraphPanel)
        hgax.title = "Visualization of the hypergraph"
        axes[hypergraphPanel] = hgax
    end

    if plot_states || plot_hyperedges || plot_active_hyperedges
        plot_count += 1
        history_box = plot_box[1, plot_count:plot_count + 1]
        plot_count += 1
        hist_plot_count = 0
    end
    
    if plot_states
        hist_plot_count += 1
        push!(panels, stateDistPanel)
        obs_lines[stateDistPanel] = []
        state_hist_box = history_box[hist_plot_count, 1]
        title = "Distribution of states"
        axes[stateDistPanel] = Axis(state_hist_box, title=title)
        num_states = length(instances(State))
        linecolors = get(colorschemes[node_colormap], 1:num_states, (1, num_states))
        for (i, state) in enumerate(instances(State))
            l = lines!(axes[stateDistPanel],
                       mo.state_history[state], 
                       label = "# of $state nodes",
                       color = linecolors[i])
            push!(obs_lines[stateDistPanel], l)
            xlims!(axes[stateDistPanel], 0, 100)
            ylims!(axes[stateDistPanel], 0, get_num_nodes(model.network))
        end
        Legend(history_box[hist_plot_count, 2],
               axes[stateDistPanel],
               framevisible = false,
               halign = :left,
               labelsize = 12)
    end
            
    if plot_hyperedges
        hist_plot_count += 1
        push!(panels, hyperedgeDistPanel)
        obs_lines[hyperedgeDistPanel] = []
        hyperedge_hist_box = history_box[hist_plot_count, 1]
        title = "Distribution of hyperdeges"
        axes[hyperedgeDistPanel] = Axis(hyperedge_hist_box, title=title)
        max_size = get_max_hyperedge_size(mo.network[])
        linecolors = get(colorschemes[hyperedge_colormap], 1:max_size, (1, max_size))
        for size in 2:max_size
            l = lines!(axes[hyperedgeDistPanel],
                       mo.hyperedge_history[size],
                       label = "hyperedges of size $size",
                       color = linecolors[size - 1])
            push!(obs_lines[hyperedgeDistPanel], l)
            xlims!(axes[hyperedgeDistPanel], 0, 100)
        end
        Legend(history_box[hist_plot_count, 2], 
               axes[hyperedgeDistPanel],
               framevisible = false,
               halign = :left,
               labelsize = 12)
    end
    
    if plot_active_hyperedges
        hist_plot_count += 1
        push!(panels, activeHyperedgesPanel)
        obs_lines[activeHyperedgesPanel] = []
        active_hist_box = history_box[hist_plot_count, 1]
        title = "Number of active hyperedges"
        axes[activeHyperedgesPanel] = Axis(active_hist_box[1, 1], title=title)
        l = lines!(axes[activeHyperedgesPanel],
                   mo.active_hyperedges_history)
        push!(obs_lines[activeHyperedgesPanel], l)
        xlims!(axes[activeHyperedgesPanel], 0, 100)
    end

    Dashboard(fig, panels, axes, obs_lines, mo, is_interactive)
end


"""
    run!(dashboard::Dashboard, num_steps::Integer, steps_per_update::Integer)

Run the simulation for `num_steps` time steps. The visualization is updated only 
once every number of steps given by `steps_per_update`.
"""
function run!(dashboard::Dashboard, num_steps::Integer, steps_per_update::Integer)
    mo = dashboard.mo
    axes = dashboard.axes

    if dashboard.is_interactive
        # TODO: something should happen here
    else
        for panel in dashboard.panels
            if panel != hypergraphPanel
                _set_hist_ax_lims!(axes[panel], num_steps)
            end
        end
        for i = 1:num_steps
            step!(mo)
            if i % steps_per_update == 0
                notify(mo.network)
                for panel in dashboard.panels
                    if panel != hypergraphPanel
                        _set_hist_ax_lims!(axes[panel], num_steps)
                    else
                        autolimits!(axes[panel])
                    end
                end
            end
        end
    end
end

function _set_hist_ax_lims!(ax::Axis, xhigh::Real)
    autolimits!(ax)
    xlims!(ax, low = 0, high=xhigh)
    ylims!(ax, low = -5)
end


"""
    record!(dashboard::Dashboard, filename::String, num_steps::Integer, steps_per_update::Integer, framerate::Integer)

Run the simulation for `num_steps` time steps and record a video of the dashboard. 

The video is saved to ./videos in a .mp4 format. `filename` should only contain the name 
of the file without the extension.
"""
function record!(dashboard::Dashboard, filename::String, num_steps::Integer, steps_per_update::Integer, framerate::Integer)
    savepath = joinpath("videos", filename * ".mp4")

    num_updates = num_steps ÷ steps_per_update
    record(dashboard.fig, savepath, 1:num_updates, framerate=framerate, compression=1) do i
        run!(dashboard, steps_per_update, steps_per_update)
    end
end


"""
    reset!(dashboard::Dashboard, model::AbstractModel)

Reset the dashboard to run the next simulation from the batch. 

The old history plot lines are made inactive and are grayed out. 
The observables in the ModelObservable are reset to track the new data from `model`.
"""
function reset!(dashboard::Dashboard, model::AbstractModel)
    mo = dashboard.mo
    axes = dashboard.axes
    # gray out the history plot lines
    if stateDistPanel in dashboard.panels
        for (i, state) in enumerate(instances(State))
            lines!(axes[stateDistPanel],
                   mo.state_history[state][],
                   linewidth = 1,
                   color = (:gray, 0.5))
        end
    end
    if hyperedgeDistPanel in dashboard.panels
        max_size = get_max_hyperedge_size(mo.network[])
        for size in 2:max_size
            lines!(axes[hyperedgeDistPanel],
                   mo.hyperedge_history[size][],
                   linewidth = 1,
                   color = (:gray, 0.5))
        end
    end
    if activeHyperedgesPanel in dashboard.panels
        lines!(axes[activeHyperedgesPanel],
               mo.active_hyperedges_history[],
               linewidth = 1,
               color = (:gray, 0.5))
    end

    # Bring the lines tied to observables in front of the gray lines
    for line in Iterators.flatten(values(dashboard.lines))
        translate!(line, 0, 0, 1)
    end

    # reset observables
    rebind_model!(mo, model)
end


"""
    save(dashboard::Dashboard, filename::String)

Saves the data from the last active run 
"""
function save(dashboard::Dashboard, filename::String)
    # TODO
end