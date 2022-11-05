export ActiveHyperedgeCount, StateCount, FinalMagnetization, HyperedgeCount, ActiveLifetime

using GLMakie

"""
A data structure for an indexed series of measurements, where both the indices and the values are Makie observables. 

This is basically a table which maps some indices (time stamps, indices, hashes, ...) to values. 
The main feature of the Log is that it can buffer values before updating the observables, since 
updating the values at every time step can be quite expensive. 
"""
mutable struct MeasurementLog{IndexType,ValueType} # TODO: rename?
    indices::Observable{Vector{IndexType}}
    values::Observable{Vector{ValueType}}
    buffered_indices::Vector{IndexType}
    buffered_values::Vector{ValueType}
    skip_points::Int64
    buffer_size::Int64
    remainder::Int64
    num_points::Int64
end

function MeasurementLog{IndexType,ValueType}(; skip_points=1,
                                             buffer_size=1) where {IndexType,ValueType}
    indices = Observable(IndexType[])
    values = Observable(ValueType[])
    buffered_indices = IndexType[]
    buffered_values = ValueType[]
    remainder = 0
    num_points = 0
    return MeasurementLog(indices, values, buffered_indices, buffered_values, skip_points,
                          buffer_size, remainder, num_points)
end

function record!(log::MeasurementLog{IndexType,ValueType}, index::IndexType,
                 value::ValueType) where {IndexType,ValueType}
    if log.buffer_size > 1
        # write the values to the buffer
        push!(log.buffered_indices, index)
        push!(log.buffered_values, value)

        # if the buffer is full, flush it 
        if length(log.buffered_indices) >= log.buffer_size
            flush_buffers!(log)
        end
    else
        # update the observables directly
        push!(log.indices[], index)
        push!(log.values[], value)
        notify(log.values)
    end
    log.num_points += 1
    return log
end

function record!(log::MeasurementLog{IndexType,ValueType},
                 value::ValueType) where {IndexType,ValueType}
    index = log.num_points
    return record!(log, index, value)
end

function flush_buffers!(log::MeasurementLog)
    # crazy mod magic to account for the fact that skip_points might not divide buffer_size
    skip = log.skip_points
    start = mod1(skip - log.remainder + 1, skip)
    log.remainder = mod(log.buffer_size - start + 1, skip)
    append!(log.indices[], log.buffered_indices[start:skip:end])
    append!(log.values[], log.buffered_values[start:skip:end])
    empty!(log.buffered_indices)
    empty!(log.buffered_values)
    notify(log.indices)
    return log
end

function save(io::IO, log::MeasurementLog)
    indices = log.indices[]
    values = log.values[]
    for (index, value) in zip(indices, values)
        println(io, "$index, $value")
    end
    return nothing
end

function clear!(log::MeasurementLog)
    empty!(log.indices[])
    empty!(log.values[])
    empty!(log.buffered_indices)
    empty!(log.buffered_values)
    log.remainder = 0
    log.num_points = 0
    return log
end

# ====================================================================================
# ------------------------------- MEASUREMENTS ----------------------------------------

abstract type AbstractMeasurement end

abstract type AbstractStepMeasurement <: AbstractMeasurement end

abstract type AbstractRunMeasurement <: AbstractMeasurement end

# ====================================================================================
# ------------------------------- Step Measurements ----------------------------------

"""
    StateCount <: AbstractStepMeasurement

Tracks the absolute number of nodes in every state.
"""
struct StateCount <: AbstractStepMeasurement
    log::MeasurementLog{Int64,Int64}
    label::State
end

function StateCount(state::State; skip_points::Int64, buffer_size::Int64)
    log = MeasurementLog{Int64,Int64}(; skip_points, buffer_size)
    return StateCount(log, state)
end

"""
    HyperedgeCount <: AbstractStepMeasurement

Tracks the number of hyperedges of every size.
"""
mutable struct HyperedgeCount <: AbstractStepMeasurement
    log::MeasurementLog{Int64,Int64}
    label::Int64
end

function HyperedgeCount(size::Int64; skip_points::Int64, buffer_size::Int64)
    log = MeasurementLog{Int64,Int64}(; skip_points, buffer_size)
    return HyperedgeCount(log, size)
end

"""
    ActiveHyperedgeCount <: AbstractStepMeasurement

Tracks the number of active hyperdeges (i.e., hyperedges with at least two nodes
in different states).
"""
mutable struct ActiveHyperedgeCount <: AbstractStepMeasurement
    log::MeasurementLog{Int64,Int64}
    label::Int64
end

function ActiveHyperedgeCount(size::Int64; skip_points::Int64, buffer_size::Int64)
    log = MeasurementLog{Int64,Int64}(; skip_points, buffer_size)
    return ActiveHyperedgeCount(log, size)
end

# ====================================================================================
# ------------------------------- Run Measurements -----------------------------------

"""
Measures the time that the system needs to deplete all active hyperedges. 
"""
struct ActiveLifetime <: AbstractRunMeasurement
    log::MeasurementLog{Int64,Int64}
end

ActiveLifetime() = ActiveLifetime(MeasurementLog{Int64,Int64}())

struct FinalMagnetization <: AbstractRunMeasurement
    log::MeasurementLog{Int64,Int64}
end

FinalMagnetization() = FinalMagnetization(MeasurementLog{Int64,Int64}())

struct FinalHyperedgeDist <: AbstractRunMeasurement
    log::MeasurementLog{Int64,Int64}
end

FinalHyperedgeDist() = FinalHyperedgeDist(MeasurementLog{Int64,Int64}())