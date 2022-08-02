function create_dashboard(model::AbstractModel,
                          history_size::Int64;
                          plot_hypergraph::Bool=false,
                          plot_states::Bool=true,
                          interactivity::Bool=false)
    f = Figure(resolution = (600, 400))
    display(f)

    mo = ModelObservable{typeof(model)}(model, history_size)

    plot_box = f[1, 1] = GridLayout()
    if interactivity
        controls_box = f[2, 1] = GridLayout()
    end

    plot_count = 0
    
    if plot_hypergraph
        plot_count += 1
        hypergraphplot(plot_box[1, plot_count], mo.network)
    end

    if plot_states
        plot_count += 1
        state_hist_box = plot_box[1, plot_count]
        for (i, state) in enumerate(instances(State))
            ax, plot = lines(state_hist_box[i, 1], mo.state_history[state])
        end
    end

    if interactivity
        # something happens
    else
        for i = 1:10
            step!(mo)
            sleep(0.1)
        end
    end
end