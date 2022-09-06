"""
    AbstractMeasurement

An abstract type for any measurements of the model.
"""
abstract type AbstractMeasurement end

"""
    AbstractTimeSeries

A time series of the evolution of some observable (as in physical observable) of the model.
"""
abstract type AbstractTimeSeries <: AbstractMeasurement end

"""
    AbstractRunMeasurement <: AbstractMeasurement

A measurement which quantifies one run of the simulation. Here, the number of 
data points is equal not to the number of time steps, but to the number of simulations 
in the batch. 
"""
abstract type AbstractRunMeasurement <: AbstractMeasurement end

# ====================================================================================
# ------------------------------- TIME SERIES ----------------------------------------

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

function record_measurement!(state_series::StateCount)
    state_dist = get_state_count(state_series.network)
    return push!(state_series.buffer, state_dist[state_series.state])
end

function record_measurement!(hyperedge_series::HyperedgeCount)
    hyperedge_dist = get_hyperedge_dist(hyperedge_series.network)
    return push!(hyperedge_series.buffer, hyperedge_dist[hyperedge_series.size])
end

function record_measurement!(active_hyperedges::ActiveHyperedgeCount)
    active_count = get_num_active_hyperedges(active_hyperedges.network)
    return push!(active_hyperedges.buffer, active_count)
end

function flush_buffers!(series::AbstractTimeSeries)
    append!(series.observable[], series.buffer)
    empty!(series.buffer)
    return series
end

function clear!(measurement::AbstractMeasurement)
    return empty!(measurement.observable[])
end

# ====================================================================================
# ------------------------------- RUN MEASUREMENTS------------------------------------

"""

Measures the time that the system needs to deplete all active hyperedges. 
"""
struct ActiveLifetime <: AbstractRunMeasurement
    observable::Observable{Vector{Int64}}
end

function ActiveLifetime()
    return ActiveLifetime(Observable(Int64[]))
end

function record_measurement!(active_lifetime::ActiveLifetime, value::Int64)
    push!(active_lifetime.observable[], value)
    notify(active_lifetime.observable)
    @show active_lifetime.observable
    return active_lifetime
end