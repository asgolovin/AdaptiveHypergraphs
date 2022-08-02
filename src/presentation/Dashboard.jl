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
end

function Dashboard(model::AbstractModel,
                   history_size::Int64;
                   plot_hypergraph::Bool=false,
                   plot_states::Bool=true,
                   interactivity::Bool=false)
    fig = Figure(resolution = (600, 400))
    display(fig)
    axes = Dict{Panel, Axis}()
    panels = Panel[]

    mo = ModelObservable{typeof(model)}(model, history_size)

    plot_box = fig[1, 1] = GridLayout()
    if interactivity
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

    if interactivity
        # something happens
    else
        for i = 1:1000
            step!(mo)
            notify(mo.network)
            if plot_states
                xlims!(axes[stateDistPanel], 0, history_size)
            end
            if plot_hypergraph
                autolimits!(axes[hypergraphPanel])
            end
            sleep(0.1)
        end
    end

    Dashboard(fig, panels, axes, mo)
end