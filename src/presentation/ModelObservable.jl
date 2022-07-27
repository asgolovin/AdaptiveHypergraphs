using DataStructures


"""
    ModelObservable{P <: PropagationRule, A <: AdaptivityRule}

Contains all observables needed to plot the evolution of the model. 
"""
struct ModelObservable{P <: PropagationRule, A <: AdaptivityRule}
    model::Observable{TS{P, A} <: AbstractModel}
    state_history::Observable{Dict{State, CircularBuffer{Int64}}}
end


function ModelObsorvable{P, A}(model::TS{P, A} <: AbstractModel) where {P <: PropagationRule,
                                                                        A <: AdaptivityRule}
    state_history = Observable(Dict{State, CircularBuffer{Int64}})
    record_state_history!(model, state_history[])
    return ModelObservable(Observale(model), state_history)
end


"""
    step!(mo::ModelObservable)
    
Progress the model one time step forward and update the history. 
"""
function step!(mo::ModelObservable)
    model, state_history = mo.model, mo.state_history
    step!(model[])
    notify(model)
    record_state_history!(model[], state_history[])
    notify(state_history)
end


"""
    record_state_history!(model::AbstractModel, state_history::Dict{State, CircularBuffer{Int64}})

Pushes the current distribution of states from the model into state_history.
"""
function record_state_history!(model::AbstractModel, 
                               state_history::Dict{State, CircularBuffer{Int64}})
    state_dist = get_state_dist(model.network)
    for state in keys(state_dist)
        push!(state_history[][state], state_dist[state])
    end
    return nothing
end
