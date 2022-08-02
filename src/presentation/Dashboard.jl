"""
Types of data that can be visualized in the dashboard.
"""
@enum Panel begin
    hypergraphPanel
    stateDistPanel
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
              is_interactive::Bool=false)

A visualization of the evolution of the hypergraph during the simulation.
"""
function Dashboard(model::AbstractModel;
                   plot_hypergraph::Bool=false,
                   plot_states::Bool=true,
                   is_interactive::Bool=false)
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
        hgax, _ = hypergraphplot(hg_box[1, 1], mo.network)
        push!(panels, hypergraphPanel)
        hgax.title = "Visualization of the hypergraph"
        axes[hypergraphPanel] = hgax
    end

    if plot_states
        plot_count += 1
        push!(panels, stateDistPanel)
        state_hist_box = plot_box[1, plot_count]
        axes[stateDistPanel] = Axis(state_hist_box[1, 1], title="Distribution of states")
        for (i, state) in enumerate(instances(State))
            lines!(axes[stateDistPanel], mo.state_history[state], label="# of $state nodes")
            ylims!(axes[stateDistPanel], 0, get_num_nodes(model.network))
        end
        state_hist_box[2, 1] = Legend(state_hist_box, 
                                      axes[stateDistPanel],
                                      orientation = :horizontal,
                                      framevisible=false)
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
                if stateDistPanel in dashboard.panels
                    xlims!(axes[stateDistPanel], 0, max(i, 100))
                end
                if hypergraphPanel in dashboard.panels
                    autolimits!(axes[hypergraphPanel])
                end
                sleep(0.1)
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