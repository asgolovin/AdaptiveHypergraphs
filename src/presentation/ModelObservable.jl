
"""
    AbstractTimeSeries

A time series of the evolution of some observable (as in physical observable) of the model.
"""
abstract type AbstractTimeSeries end

"""
    StateCount <: AbstractTimeSeries

Tracks the absolute number of nodes in state `state`. 
"""
mutable struct StateCount <: AbstractTimeSeries
    network::HyperNetwork
    state::State
    observable::Observable{Vector{Int64}}
    buffer::Vector{Int64}
end

function StateCount(network::HyperNetwork, state::State)
    return StateCount(network, state, Observable(Int64[]), Int64[])
end

"""
    HyperedgeCount <: AbstractTimeSeries

Tracks the number of hyperedges of size `size`. 
"""
mutable struct HyperedgeCount <: AbstractTimeSeries
    network::HyperNetwork
    size::Int64
    observable::Observable{Vector{Int64}}
    buffer::Vector{Int64}
end

function HyperedgeCount(network::HyperNetwork, size::Int64)
    return HyperedgeCount(network, size, Observable(Int64[]), Int64[])
end

"""
    ActiveHyperedgeCount <: AbstractTimeSeries

Tracks the number of active hyperdeges (i.e., hyperedges with at least two nodes
in different states).
"""
mutable struct ActiveHyperedgeCount <: AbstractTimeSeries
    network::HyperNetwork
    observable::Observable{Vector{Int64}}
    buffer::Vector{Int64}
end

function ActiveHyperedgeCount(network::HyperNetwork)
    return ActiveHyperedgeCount(network, Observable(Int64[]), Int64[])
end

function record_history!(state_series::StateCount)
    state_dist = get_state_count(state_series.network)
    return push!(state_series.buffer, state_dist[state_series.state])
end

function record_history!(hyperedge_series::HyperedgeCount)
    hyperedge_dist = get_hyperedge_dist(hyperedge_series.network)
    return push!(hyperedge_series.buffer, hyperedge_dist[hyperedge_series.size])
end

function record_history!(active_hyperedges::ActiveHyperedgeCount)
    active_count = get_num_active_hyperedges(active_hyperedges.network)
    return push!(active_hyperedges.buffer, active_count)
end

function flush_buffers!(series::AbstractTimeSeries)
    append!(series.observable[], series.buffer)
    empty!(series.buffer)
    notify(series.observable)
    return series
end

function clear!(series::AbstractTimeSeries)
    return empty!(series.observable[])
end

"""
    ModelObservable{M <: AbstarctModel}

Contains all observables (both as in physical observables and as in Makie Observable class) 
needed to plot the evolution of the model. 

This class be methaphorically thought of as a scientist who runs some experiments 
(i.e., pushes a button to evolve the model one step forward), observes the results 
and puts a new marker into a plot. ModelObservable provides its own interface to the 
step! function which evolves the model and records the results. 
"""
struct ModelObservable{M<:AbstractModel}
    model::Observable{M}
    network::Observable{HyperNetwork}
    state_series::Vector{StateCount}
    hyperedge_series::Vector{HyperedgeCount}
    active_hyperedges_series::ActiveHyperedgeCount

    function ModelObservable{M}(model::M) where {M<:AbstractModel}
        state_series = [StateCount(model.network, state) for state in instances(State)]
        max_size = get_max_hyperedge_size(model.network)
        hyperedge_series = [HyperedgeCount(model.network, size) for size in 2:max_size]
        active_hyperedges_series = ActiveHyperedgeCount(model.network)

        mo = new{M}(Observable(model),
                    Observable(model.network),
                    state_series,
                    hyperedge_series,
                    active_hyperedges_series)
        return record_history!(mo)
    end
end

function _get_all_series(mo::ModelObservable)
    return [mo.state_series; mo.hyperedge_series; mo.active_hyperedges_series]
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
    for series in _get_all_series(mo)
        record_history!(series)
    end
    return mo
end

function flush_buffers!(mo::ModelObservable)
    for series in _get_all_series(mo)
        flush_buffers!(series)
    end
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
    for series in _get_all_series(mo)
        series.network = model.network
    end
    return record_history!(mo)
end

"""
    clear!(mo::ModelObservable)

Clear the history buffers in the ModelObservable.
"""
function clear!(mo::ModelObservable)
    for series in _get_all_series(mo)
        clear!(series)
    end
    return mo
end