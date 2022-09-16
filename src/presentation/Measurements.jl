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
    time_steps::Observable{Vector{Int64}}
    values::Observable{Vector{Int64}}
    buffer::Vector{Int64}
    skip_points::Int64
end

function StateCount(network::HyperNetwork, state::State, skip_points::Int64)
    return StateCount(network, state, Observable(Int64[]), Observable(Int64[]),
                      Int64[], skip_points)
end

"""
    HyperedgeCount <: AbstractTimeSeries

Tracks the number of hyperedges of size `size`. 
"""
mutable struct HyperedgeCount <: AbstractTimeSeries
    network::HyperNetwork
    size::Int64
    time_steps::Observable{Vector{Int64}}
    values::Observable{Vector{Int64}}
    buffer::Vector{Int64}
    skip_points::Int64
end

function HyperedgeCount(network::HyperNetwork, size::Int64, skip_points::Int64)
    return HyperedgeCount(network, size, Observable(Int64[]), Observable(Int64[]),
                          Int64[], skip_points)
end

"""
    ActiveHyperedgeCount <: AbstractTimeSeries

Tracks the number of active hyperdeges (i.e., hyperedges with at least two nodes
in different states).
"""
mutable struct ActiveHyperedgeCount <: AbstractTimeSeries
    network::HyperNetwork
    size::Int64
    time_steps::Observable{Vector{Int64}}
    values::Observable{Vector{Int64}}
    buffer::Vector{Int64}
    skip_points::Int64
end

function ActiveHyperedgeCount(network::HyperNetwork, size::Int64, skip_points::Int64)
    return ActiveHyperedgeCount(network, size, Observable(Int64[]), Observable(Int64[]),
                                Int64[], skip_points)
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
    active_count = get_num_active_hyperedges(active_hyperedges.network,
                                             active_hyperedges.size)
    return push!(active_hyperedges.buffer, active_count)
end

function flush_buffers!(series::AbstractTimeSeries)
    skip = series.skip_points
    append!(series.values[], series.buffer[1:skip:end])
    start_time = length(series.time_steps[]) == 0 ? 1 : series.time_steps[][end] + skip
    end_time = start_time + length(series.buffer) - 1
    append!(series.time_steps[], start_time:skip:end_time)
    empty!(series.buffer)
    return series
end

function clear!(measurement::AbstractMeasurement)
    empty!(measurement.values[])
    empty!(measurement.time_steps[])
    return measurement
end

# ====================================================================================
# ------------------------------- RUN MEASUREMENTS------------------------------------

"""

Measures the time that the system needs to deplete all active hyperedges. 
"""
struct ActiveLifetime <: AbstractRunMeasurement
    values::Observable{Vector{Int64}}
end

function ActiveLifetime()
    return ActiveLifetime(Observable(Int64[]))
end

function record_measurement!(active_lifetime::ActiveLifetime, value::Int64)
    push!(active_lifetime.values[], value)
    notify(active_lifetime.values)
    return active_lifetime
end

struct FinalMagnetization <: AbstractRunMeasurement
    values::Observable{Vector{Int64}}
    has_converged::Observable{Vector{Bool}}
end

function FinalMagnetization()
    return FinalMagnetization(Observable(Int64[]), Observable(Bool[]))
end

function record_measurement!(final_magnetization::FinalMagnetization, value::Int64,
                             has_converged::Bool)
    push!(final_magnetization.values[], value)
    push!(final_magnetization.has_converged[], has_converged)
    notify(final_magnetization.values)
    return final_magnetization
end