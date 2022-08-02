using DataStructures


"""
    ModelObservable{P <: PropagationRule, A <: AdaptivityRule}

Contains all observables needed to plot the evolution of the model. 
"""
struct ModelObservable{M <: AbstractModel}
    model::Observable{M}
    network::Observable{HyperNetwork}
    state_history::Dict{State, Observable{Vector{Int64}}}
end


function ModelObservable{M}(model::M) where {M <: AbstractModel}
    state_history = Dict{State, Observable{Vector{Int64}}}()
    for state in instances(State)
        state_history[state] = Observable(Vector{Int64}())
    end
    record_state_history!(model, state_history)
    return ModelObservable(Observable(model), 
                           Observable(model.network),
                           state_history)
end


"""
    step!(mo::ModelObservable)
    
Progress the model one time step forward and update the history. 
"""
function step!(mo::ModelObservable)
    model, state_history = mo.model, mo.state_history
    network_changed = step!(model[])
    network_changed && notify(model)
    record_state_history!(model[], state_history)
end


"""
    record_state_history!(model::AbstractModel, state_history::Dict{State, Vector{Int64}})

Pushes the current distribution of states from the model into state_history.
"""
function record_state_history!(model::AbstractModel, 
                               state_history::Dict{State, Observable{Vector{Int64}}})
    state_dist = get_state_dist(model.network)
    for state in keys(state_dist)
        push!(state_history[state][], state_dist[state])
        notify(state_history[state])
    end
    return nothing
end
