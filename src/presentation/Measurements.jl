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
    return series
end

function clear!(series::AbstractTimeSeries)
    return empty!(series.observable[])
end
