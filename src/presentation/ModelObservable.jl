"""
    ModelObservable{P <: PropagationRule, A <: AdaptivityRule}

Contains all observables needed to plot the evolution of the model. 
"""
struct ModelObservable{M<:AbstractModel}
    model::Observable{M}
    network::Observable{HyperNetwork}
    state_history::Dict{State,Observable{Vector{Int64}}}
    hyperedge_history::Dict{Int64,Observable{Vector{Int64}}}
    active_hyperedges_history::Observable{Vector{Int64}}

    function ModelObservable{M}(model::M) where {M<:AbstractModel}
        state_history = Dict{State,Observable{Vector{Int64}}}()
        hyperedge_history = Dict{Int64,Observable{Vector{Int64}}}()
        active_hyperedges_history = Observable(Vector{Int64}())
        for state in instances(State)
            state_history[state] = Observable(Vector{Int64}())
        end
        for size in 2:get_max_hyperedge_size(model.network)
            hyperedge_history[size] = Observable(Vector{Int64}())
        end
        mo = new{M}(Observable(model),
                    Observable(model.network),
                    state_history,
                    hyperedge_history,
                    active_hyperedges_history)
        return record_history!(mo)
    end
end

"""
    step!(mo::ModelObservable)
    
Progress the model one time step forward and update the history. 
"""
function step!(mo::ModelObservable)
    network_changed = step!(mo.model[])
    network_changed && notify(mo.model)
    return record_history!(mo)
end

"""
    record_history!(mo::ModelObservable)

Push the current distribution of states from the model into state_history.
"""
function record_history!(mo::ModelObservable)
    state_history = mo.state_history
    hyperedge_history = mo.hyperedge_history
    active_hyperedges_history = mo.active_hyperedges_history
    network = mo.network[]
    state_dist = get_state_count(network)

    for state in keys(state_dist)
        push!(state_history[state][], state_dist[state])
        notify(state_history[state])
    end
    hyperedge_dist = get_hyperedge_dist(network)
    for size in keys(hyperedge_dist)
        push!(hyperedge_history[size][], hyperedge_dist[size])
        notify(hyperedge_history[size])
    end
    active_count = 0
    for hyperedge in get_hyperedges(network)
        if is_active(network, hyperedge)
            active_count += 1
        end
    end
    push!(active_hyperedges_history[], active_count)
    notify(active_hyperedges_history)
    return mo
end

"""
    rebind_model!(mo::ModelObservable, model::AbstractModel)

Rebind the observables to track the new `model`.
"""
function rebind_model!(mo::ModelObservable, model::AbstractModel)
    clear!(mo)
    mo.model[] = model
    mo.network[] = model.network
    return record_history!(mo)
end

"""
    clear!(mo::ModelObservable)

Clear the history buffers in the ModelObservable.
"""
function clear!(mo::ModelObservable)
    mo.active_hyperedges_history[] = Vector{Int64}()
    for state in instances(State)
        mo.state_history[state][] = Vector{Int64}()
    end
    for size in 2:get_max_hyperedge_size(mo.network[])
        mo.hyperedge_history[size][] = Vector{Int64}()
    end
    return mo
end