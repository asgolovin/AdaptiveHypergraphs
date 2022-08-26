"""
    ModelObservable{P <: PropagationRule, A <: AdaptivityRule}

Contains all observables needed to plot the evolution of the model. 
"""
struct ModelObservable{M<:AbstractModel}
    model::Observable{M}
    network::Observable{HyperNetwork}
    state_history::Dict{State,Observable{Vector{Int64}}}
    state_history_buffer::Dict{State,Vector{Int64}}
    hyperedge_history::Dict{Int64,Observable{Vector{Int64}}}
    hyperedge_history_buffer::Dict{Int64,Vector{Int64}}
    active_hyperedges_history::Observable{Vector{Int64}}
    active_hyperedges_history_buffer::Vector{Int64}

    function ModelObservable{M}(model::M) where {M<:AbstractModel}
        state_history = Dict{State,Observable{Vector{Int64}}}()
        state_history_buffer = Dict{State,Vector{Int64}}()
        hyperedge_history = Dict{Int64,Observable{Vector{Int64}}}()
        hyperedge_history_buffer = Dict{Int64,Vector{Int64}}()
        active_hyperedges_history = Observable(Vector{Int64}())
        active_hyperedges_history_buffer = Int64[]
        for state in instances(State)
            state_history[state] = Observable(Vector{Int64}())
            state_history_buffer[state] = Int64[]
        end
        for size in 2:get_max_hyperedge_size(model.network)
            hyperedge_history[size] = Observable(Vector{Int64}())
            hyperedge_history_buffer[size] = Int64[]
        end
        mo = new{M}(Observable(model),
                    Observable(model.network),
                    state_history,
                    state_history_buffer,
                    hyperedge_history,
                    hyperedge_history_buffer,
                    active_hyperedges_history,
                    active_hyperedges_history_buffer)
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

Push the current distribution of states from the model into the history buffer vectors. 
"""
function record_history!(mo::ModelObservable)
    state_dist = get_state_count(mo.network[])
    for state in keys(state_dist)
        push!(mo.state_history_buffer[state], state_dist[state])
    end

    hyperedge_dist = get_hyperedge_dist(mo.network[])
    for size in keys(hyperedge_dist)
        push!(mo.hyperedge_history_buffer[size], hyperedge_dist[size])
    end
    active_count = get_num_active_hyperedges(mo.network[])
    push!(mo.active_hyperedges_history_buffer, active_count)
    return mo
end

function flush_buffers!(mo::ModelObservable)
    for state in keys(mo.state_history_buffer)
        append!(mo.state_history[state][], mo.state_history_buffer[state])
        empty!(mo.state_history_buffer[state])
        notify(mo.state_history[state])
    end
    for size in keys(mo.hyperedge_history_buffer)
        append!(mo.hyperedge_history[size][], mo.hyperedge_history_buffer[size])
        empty!(mo.hyperedge_history_buffer[size])
        notify(mo.hyperedge_history[size])
    end
    append!(mo.active_hyperedges_history[], mo.active_hyperedges_history_buffer)
    empty!(mo.active_hyperedges_history_buffer)
    notify(mo.active_hyperedges_history)
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