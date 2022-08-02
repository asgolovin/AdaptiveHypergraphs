"""
    ModelObservable{P <: PropagationRule, A <: AdaptivityRule}

Contains all observables needed to plot the evolution of the model. 
"""
struct ModelObservable{M <: AbstractModel}
    model::Observable{M}
    network::Observable{HyperNetwork}
    state_history::Dict{State, Observable{Vector{Int64}}}
    hyperedge_history::Dict{Int64, Observable{Vector{Int64}}}
end


function ModelObservable{M}(model::M) where {M <: AbstractModel}
    state_history = Dict{State, Observable{Vector{Int64}}}()
    hyperedge_history = Dict{Int64, Observable{Vector{Int64}}}()
    for state in instances(State)
        state_history[state] = Observable(Vector{Int64}())
    end
    for size in 2:get_max_hyperedge_size(model.network)
        hyperedge_history[size] = Observable(Vector{Int64}())
    end
    record_history!(model, state_history, hyperedge_history)
    return ModelObservable(Observable(model), 
                           Observable(model.network),
                           state_history,
                           hyperedge_history)
end


"""
    step!(mo::ModelObservable)
    
Progress the model one time step forward and update the history. 
"""
function step!(mo::ModelObservable)
    network_changed = step!(mo.model[])
    network_changed && notify(mo.model)
    record_history!(mo.model[], mo.state_history, mo.hyperedge_history)
end


"""
    record_history!(model::AbstractModel, state_history::Dict{State, Observable{Vector{Int64}}}, hyperedge_history::Dict{Int64, Observable{Vector{Int64}}})

Pushes the current distribution of states from the model into state_history.
"""
function record_history!(model::AbstractModel, 
                         state_history::Dict{State, Observable{Vector{Int64}}},
                         hyperedge_history::Dict{Int64, Observable{Vector{Int64}}})
    state_dist = get_state_dist(model.network)
    for state in keys(state_dist)
        push!(state_history[state][], state_dist[state])
        notify(state_history[state])
    end
    hyperedge_dist = get_hyperedge_dist(model.network)
    for size in 2:get_max_hyperedge_size(model.network)
        push!(hyperedge_history[size][], hyperedge_dist[size])
        notify(hyperedge_history[size])
    end
    return nothing
end
