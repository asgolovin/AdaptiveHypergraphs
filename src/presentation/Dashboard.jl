"""
Types of data that can be visualized in the dashboard.

Each panel can have multiple axes: for example, the state_dist_panel 
"""
@enum Panel begin
    hypergraphPanel
    stateDistPanel
end

struct Dashboard
    fig::Figure
    panels::Vector{Panel}
    axes::Dict{Panel, Vector{Axis}}
    mo::ModelObservable
end

function Dashboard(model::AbstractModel,
                   history_size::Int64;
                   plot_hypergraph::Bool=false,
                   plot_states::Bool=true,
                   interactivity::Bool=false)
    fig = Figure(resolution = (600, 400))
    display(fig)
    axes = Dict{Panel, Vector{Axis}}()
    panels = Panel[]

    mo = ModelObservable{typeof(model)}(model, history_size)

    plot_box = fig[1, 1] = GridLayout()
    if interactivity
        controls_box = fig[2, 1] = GridLayout()
    end

    plot_count = 0
    
    if plot_hypergraph
        plot_count += 1
        hgax, _ = hypergraphplot(plot_box[1, plot_count], mo.network)
        push!(panels, hypergraphPanel)
        axes[hypergraphPanel] = [hgax, ]
    end

    if plot_states
        plot_count += 1
        push!(panels, stateDistPanel)
        axes[stateDistPanel] = []
        state_hist_box = plot_box[1, plot_count]
        for (i, state) in enumerate(instances(State))
            ax, _ = lines(state_hist_box[i, 1], mo.state_history[state])
            push!(axes[stateDistPanel], ax)
        end

    end

    if interactivity
        # something happens
    else
        for i = 1:100
            step!(mo)
            notify(mo.network)
            if plot_states
                for (i, state) in enumerate(instances(State))
                    autolimits!(axes[stateDistPanel][i])
                end
            end
            sleep(0.1)
        end
    end

    Dashboard(fig, panels, axes, mo)
end