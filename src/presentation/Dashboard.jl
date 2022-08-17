export Panel, Dashboard, run!, record!

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
    mo::ModelObservable
    is_interactive::Bool
end


"""
    Dashboard(model::AbstractModel;
              plot_hypergraph::Bool=false,
              plot_states::Bool=true,
              plot_hyperedges::Bool=true,
              plot_active_hyperedges::Bool=true,
              is_interactive::Bool=false)

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

    if plot_states || plot_hyperedges
        plot_count += 1
        history_box = plot_box[1, plot_count]
    end

    if plot_states
        push!(panels, stateDistPanel)
        state_hist_box = history_box[1, 1]
        axes[stateDistPanel] = Axis(state_hist_box[1, 1], title="Distribution of states")
        num_states = length(instances(State))
        linecolors = get(colorschemes[node_colormap], 1:num_states, (1, num_states))
        for (i, state) in enumerate(instances(State))
            lines!(axes[stateDistPanel],
                   mo.state_history[state], 
                   label = "# of $state nodes",
                   color = linecolors[i])
            xlims!(axes[stateDistPanel], 0, 100)
            ylims!(axes[stateDistPanel], 0, get_num_nodes(model.network))
        end
        state_hist_box[2, 1] = Legend(state_hist_box, 
                                      axes[stateDistPanel],
                                      orientation = :horizontal,
                                      framevisible=false)
    end

    if plot_hyperedges
        push!(panels, hyperedgeDistPanel)
        hyperedge_hist_box = history_box[plot_states ? 2 : 1, 1]
        axes[hyperedgeDistPanel] = Axis(hyperedge_hist_box[1, 1], title="Distribution of hyperdeges")
        max_hyperedge_size = get_max_hyperedge_size(mo.network[])
        linecolors = get(colorschemes[hyperedge_colormap], 1:max_hyperedge_size, (1, max_hyperedge_size))
        for size in 2:max_hyperedge_size
            lines!(axes[hyperedgeDistPanel],
                   mo.hyperedge_history[size],
                   label="# of hyperdeges of size $size",
                   color = linecolors[size - 1])
            xlims!(axes[hyperedgeDistPanel], 0, 100)
        end
        hyperedge_hist_box[2, 1] = Legend(hyperedge_hist_box, 
                                          axes[hyperedgeDistPanel],
                                          orientation = :vertical,
                                          framevisible=false)
    end

    if plot_active_hyperedges
        push!(panels, activeHyperedgesPanel)
        plot_count += 1
        active_hist_box = plot_box[1, plot_count]
        axes[activeHyperedgesPanel] = Axis(active_hist_box[1, 1], title="Number of active hyperedges")
        lines!(axes[activeHyperedgesPanel],
               mo.active_hyperedges_history)
        xlims!(axes[activeHyperedgesPanel], 0, 100)
    end

    Dashboard(fig, panels, axes, mo, is_interactive)
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
        for i = 1:num_steps
            step!(mo)
            if i % steps_per_update == 0
                notify(mo.network)
                for panel in dashboard.panels
                    if panel != hypergraphPanel
                        autolimits!(axes[panel])
                        xlims!(axes[panel], 0, max(i, 100))
                        ylims!(axes[panel], low = 0)
                    else
                        autolimits!(axes[panel])
                    end
                end
                sleep(0.5)
            end
        end
    end
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