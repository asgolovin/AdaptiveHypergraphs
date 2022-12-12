export ActiveHyperedgeCount, StateCount, FinalMagnetization, HyperedgeCount, ActiveLifetime,
       AvgHyperedgeCount, SlowManifoldFit

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
    auto_notify::Bool
    remainder::Int64
    num_points::Int64
end

"""
    function MeasurementLog{IndexType,ValueType}(; skip_points=1,
        buffer_size=1,
        auto_notify=true) where {IndexType,ValueType}

# Arguments
- `skip_points::Integer = 1` - if this is greater than one, then only every n-th point gets written to the observables. Only works if `buffer_size` > 1
- `buffer_size::Integer = 0` - if this is greater than zero, then the points are written to a non-Observable buffer of given size. Once the buffer is full, they are pushed to the actual Observable. This improves performance if the updates are very frequent.
- `auto_notify::Bool = true` - if set to true, any operations which modify the observables will automatically notify them. In the other case, this has to be done explicitly using notify(log). This can be used in cases where plots depend on multiple different logs and it is important to ensure that all logs are updated with the new values before the observables are triggered. 
"""
function MeasurementLog{IndexType,ValueType}(; skip_points=1,
                                             buffer_size=0,
                                             auto_notify=true) where {IndexType,ValueType}
    indices = Observable(IndexType[])
    values = Observable(ValueType[])
    buffered_indices = IndexType[]
    buffered_values = ValueType[]
    remainder = 0
    num_points = 0
    return MeasurementLog(indices, values, buffered_indices, buffered_values, skip_points,
                          buffer_size, auto_notify, remainder, num_points)
end

function MeasurementLog{IndexType,ValueType}(indices::Observable{Vector{IndexType}},
                                             values::Observable{Vector{ValueType}}) where {IndexType,
                                                                                           ValueType}
    @assert length(indices[]) == length(values[])
    buffered_indices = IndexType[]
    buffered_values = ValueType[]
    skip_points = 1
    buffer_size = 0
    auto_notify = true
    remainder = 0
    num_points = length(indices[])
    return MeasurementLog(indices, values, buffered_indices, buffered_values, skip_points,
                          buffer_size, auto_notify, remainder, num_points)
end

function Base.show(io::IO, log::MeasurementLog)
    return print(io, "$(typeof(log))()")
end

function Base.show(io::IO, ::MIME"text/plain", log::MeasurementLog)
    println(io, log)

    indices = log.indices[]
    values = log.values[]
    if length(log.indices[]) > 10
        first_indices = join(indices[1:5], ", ")
        first_values = join(values[1:5], ", ")
        last_indices = join(indices[(end - 4):end], ", ")
        last_values = join(values[(end - 4):end], ", ")
        indices = "[$first_indices, ..., $last_indices]"
        values = "[$first_values, ..., $last_values]"
    end
    println(io, "  indices: $indices")
    println(io, "  values: $values")
    return println(io, "  buffer: $(length(log.buffered_indices))/$(log.buffer_size) full")
end

"""
    record!(log::MeasurementLog{IndexType,ValueType}, index::IndexType,
        value::ValueType) where {IndexType,ValueType}

Push a new pair of (`index`, `value`) into the log. 
"""
function record!(log::MeasurementLog{IndexType,ValueType}, index::IndexType,
                 value::ValueType) where {IndexType,ValueType}
    if log.buffer_size > 0
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
        if log.auto_notify
            notify(log.indices)
            notify(log.values)
        end
    end
    log.num_points += 1
    return log
end

function record!(log::MeasurementLog{IndexType,ValueType},
                 value::ValueType) where {IndexType,ValueType}
    index = log.num_points
    return record!(log, index, value)
end

"""
    flush_buffers!(log::MeasurementLog)

Write the values in the buffers to the observables. If `auto_notify` is set to `true`, the observables are notified after the update. 
"""
function flush_buffers!(log::MeasurementLog)
    # crazy mod magic to account for the fact that skip_points might not divide buffer_size
    skip = log.skip_points
    start = mod1(skip - log.remainder + 1, skip)
    log.remainder = mod(log.buffer_size - start + 1, skip)
    append!(log.indices[], log.buffered_indices[start:skip:end])
    append!(log.values[], log.buffered_values[start:skip:end])
    empty!(log.buffered_indices)
    empty!(log.buffered_values)
    if log.auto_notify
        notify(log.indices)
        notify(log.values)
    end
    return log
end

function GLMakie.notify(log::MeasurementLog)
    GLMakie.notify(log.indices)
    GLMakie.notify(log.values)
    return nothing
end

"""
    save(io::IO, log::MeasurementLog)

Write the log into the IO stream as a CSV table separated by comma. 
"""
function save(io::IO, log::MeasurementLog)
    indices = log.indices[]
    values = log.values[]
    for (index, value) in zip(indices, values)
        println(io, "$index, $value")
    end
    return nothing
end

"""
    clear(log::MeasurementLog)

Empty all vectors and reset the log. 
"""
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

abstract type AbstractBatchMeasurement <: AbstractMeasurement end

function Base.getproperty(obj::AbstractMeasurement, sym::Symbol)
    # short-hand to access the log information. 
    if sym === :values
        return obj.log.values
    elseif sym === :indices
        return obj.log.indices
    else
        return getfield(obj, sym)
    end
end

function Base.show(io::IO, meas::AbstractMeasurement)
    if :label in propertynames(meas)
        arguments = "$(meas.label)"
    else
        arguments = ""
    end
    return print(io, "$(typeof(meas))($arguments)")
end

function Base.show(io::IO, ::MIME"text/plain", meas::AbstractMeasurement)
    println(io, meas)
    print(io, " log: ")
    return show(io, MIME("text/plain"), meas.log)
end

# ====================================================================================
# ------------------------------- Step Measurements ----------------------------------

"""
    StateCount <: AbstractStepMeasurement

Tracks the absolute number of nodes in every state.
"""
struct StateCount <: AbstractStepMeasurement
    log::MeasurementLog{Float64,Int64}
    label::State
end

function StateCount(state::State; skip_points::Int64=1, buffer_size::Int64=0)
    log = MeasurementLog{Float64,Int64}(; skip_points, buffer_size, auto_notify=false)
    return StateCount(log, state)
end

"""
    HyperedgeCount <: AbstractStepMeasurement

Tracks the number of hyperedges of every size.
"""
mutable struct HyperedgeCount <: AbstractStepMeasurement
    log::MeasurementLog{Float64,Int64}
    label::Int64
end

function HyperedgeCount(size::Int64; skip_points::Int64=1, buffer_size::Int64=0)
    log = MeasurementLog{Float64,Int64}(; skip_points, buffer_size, auto_notify=false)
    return HyperedgeCount(log, size)
end

"""
    ActiveHyperedgeCount <: AbstractStepMeasurement

Tracks the number of active hyperdeges (i.e., hyperedges with at least two nodes
in different states).
"""
mutable struct ActiveHyperedgeCount <: AbstractStepMeasurement
    log::MeasurementLog{Float64,Int64}
    label::Int64
end

function ActiveHyperedgeCount(size::Int64; skip_points::Int64=1, buffer_size::Int64=0)
    log = MeasurementLog{Float64,Int64}(; skip_points, buffer_size, auto_notify=false)
    return ActiveHyperedgeCount(log, size)
end

# ====================================================================================
# ------------------------------- Run Measurements -----------------------------------

"""
Measures the time that the system needs to deplete all active hyperedges. 
"""
struct ActiveLifetime <: AbstractRunMeasurement
    log::MeasurementLog{Int64,Float64}
end

ActiveLifetime() = ActiveLifetime(MeasurementLog{Int64,Float64}())

"""
Measures the magnetization in the end of the simulation. 
"""
struct FinalMagnetization <: AbstractRunMeasurement
    log::MeasurementLog{Int64,Int64}
end

FinalMagnetization() = FinalMagnetization(MeasurementLog{Int64,Int64}())

struct AvgHyperedgeCount <: AbstractRunMeasurement
    log::MeasurementLog{Int64,Float64}
    label::Int64
end

function AvgHyperedgeCount(size::Int64)
    return AvgHyperedgeCount(MeasurementLog{Int64,Float64}(), size)
end

struct SlowManifoldFit <: AbstractRunMeasurement
    log::MeasurementLog{Int64,NTuple{3,Float64}}
    label::Int64
end

function SlowManifoldFit(size::Int64)
    log = MeasurementLog{Int64,NTuple{3,Float64}}()
    return SlowManifoldFit(log, size)
end