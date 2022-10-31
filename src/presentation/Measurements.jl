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
            # crazy mod magic to account for the fact that skip_points might not divide buffer_size
            skip = log.skip_points
            start = mod1(skip - log.remainder + 1, skip)
            log.remainder = mod(log.buffer_size - start + 1, skip)
            append!(log.indices[], log.buffered_indices[start:skip:end])
            append!(log.values[], log.buffered_values[start:skip:end])
            empty!(log.buffered_indices)
            empty!(log.buffered_values)
            notify(log.values)
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
    return log
end

# ====================================================================================
# ------------------------------- MEASUREMENTS ----------------------------------------

abstract type AbstractMeasurement end

abstract type AbstractStepMeasurement end

abstract type AbstractRunMeasurement end

# ====================================================================================
# ------------------------------- TIME SERIES ----------------------------------------

"""
    StateCount <: AbstractStepMeasurement

Tracks the absolute number of nodes in state `state`. 
"""
struct StateCount <: AbstractStepMeasurement
    log::Dict{State,MeasurementLog{Int64,Int64}}
end

function StateCount(skip_points::Int64, buffer_size::Int64)
    log = Dict(state => MeasurementLog{Int64,Int64}(; skip_points, buffer_size)
               for state in instances(State))
    return StateCount(log)
end

"""
    HyperedgeCount <: AbstractStepMeasurement

Tracks the number of hyperedges of size `size`. 
"""
mutable struct HyperedgeCount <: AbstractStepMeasurement
    log::Dict{Int64,MeasurementLog{Int64,Int64}}
end

function HyperedgeCount(max_size::Int64, skip_points::Int64, buffer_size::Int64)
    log = Dict(size => MeasurementLog{Int64,Int64}(; skip_points, buffer_size)
               for size in 2:max_size)
    return HyperedgeCount(log)
end

"""
    ActiveHyperedgeCount <: AbstractStepMeasurement

Tracks the number of active hyperdeges (i.e., hyperedges with at least two nodes
in different states).
"""
mutable struct ActiveHyperedgeCount <: AbstractStepMeasurement
    log::Dict{Int64,MeasurementLog{Int64,Int64}}
end

function ActiveHyperedgeCount(max_size::Int64, skip_points::Int64, buffer_size::Int64)
    log = Dict(size => MeasurementLog{Int64,Int64}(; skip_points, buffer_size)
               for size in 2:max_size)
    return ActiveHyperedgeCount(log)
end

# ====================================================================================
# ------------------------------- RUN MEASUREMENTS------------------------------------

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