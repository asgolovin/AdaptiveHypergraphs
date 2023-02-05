export MeasurementLog, ActiveHyperedgeCount, StateCount, MotifCount, FinalMagnetization,
       HyperedgeCount,
       ActiveLifetime, FakeDiffEq, AvgHyperedgeCount, SlowManifoldFit

using Observables

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
    last_value::Union{ValueType,Nothing}
    second_to_last_value::Union{ValueType,Nothing}
    save_file::Union{Nothing,String}
end

"""
    function MeasurementLog{IndexType,ValueType}(; skip_points=1,
        buffer_size=1,
        auto_notify=true) where {IndexType,ValueType}

# Arguments
- `skip_points::Int64 = 1` - if this is greater than one, then only every n-th point gets written to the observables. Only works if `buffer_size` > 1
- `buffer_size::Int64 = 0` - if this is greater than zero, then the points are written to a non-Observable buffer of given size. Once the buffer is full, they are pushed to the actual Observable. This improves performance if the updates are very frequent.
- `auto_notify::Bool = true` - if set to true, any operations which modify the observables will automatically notify them. In the other case, this has to be done explicitly using notify(log). This can be used in cases where plots depend on multiple different logs and it is important to ensure that all logs are updated with the new values before the observables are triggered. 
"""
function MeasurementLog{IndexType,ValueType}(; skip_points::Int64=1,
                                             buffer_size::Int64=0,
                                             save_file::Union{Nothing,String}=nothing,
                                             auto_notify::Bool=true) where {IndexType,
                                                                            ValueType}
    indices = Observable(IndexType[])
    values = Observable(ValueType[])
    buffered_indices = IndexType[]
    buffered_values = ValueType[]
    remainder = 0
    num_points = 0
    last_value = nothing
    second_to_last_value = nothing
    return MeasurementLog(indices, values, buffered_indices, buffered_values, skip_points,
                          buffer_size, auto_notify, remainder, num_points, last_value,
                          second_to_last_value, save_file)
end

function MeasurementLog{IndexType,ValueType}(indices::Observable{Vector{IndexType}},
                                             values::Observable{Vector{ValueType}};
                                             save_file::Union{Nothing,String}=nothing) where
         {IndexType,
          ValueType}
    @assert length(indices[]) == length(values[])
    buffered_indices = IndexType[]
    buffered_values = ValueType[]
    skip_points = 1
    buffer_size = 0
    auto_notify = true
    remainder = 0
    num_points = length(indices[])
    last_value = nothing
    second_to_last_value = nothing
    return MeasurementLog(indices, values, buffered_indices, buffered_values, skip_points,
                          buffer_size, auto_notify, remainder, num_points, last_value,
                          second_to_last_value, save_file)
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
    log.second_to_last_value = log.last_value
    log.last_value = value
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
        # write to file
        if typeof(log.save_file) <: String
            open(log.save_file, "a") do io
                return write(io, "$(log.num_points + 1), $index, $value\n")
            end
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
    # write the contents of the buffer to file
    if typeof(log.save_file) <: String
        open(log.save_file, "a") do io
            point_id = log.num_points - length(log.buffered_indices) + 1
            for (index, value) in zip(log.buffered_indices, log.buffered_values)
                write(io, "$point_id, $index, $value\n")
                point_id += 1
            end
        end
    end

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

function Observables.notify(log::MeasurementLog)
    Observables.notify(log.indices)
    Observables.notify(log.values)
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

function StateCount(state::State; skip_points::Int64=1, buffer_size::Int64=0,
                    save_folder::Union{Nothing,String}=nothing)
    save_file = _create_save_file(save_folder, "state_count", state)
    log = MeasurementLog{Float64,Int64}(; skip_points, buffer_size, save_file,
                                        auto_notify=false)
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

function HyperedgeCount(size::Int64; skip_points::Int64=1, buffer_size::Int64=0,
                        save_folder::Union{Nothing,String})
    save_file = _create_save_file(save_folder, "hyperedge_count", size)
    log = MeasurementLog{Float64,Int64}(; skip_points, buffer_size, save_file,
                                        auto_notify=false)
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

function ActiveHyperedgeCount(size::Int64; skip_points::Int64=1, buffer_size::Int64=0,
                              save_folder::Union{Nothing,String})
    save_file = _create_save_file(save_folder, "active_hyperedge_count", size)
    log = MeasurementLog{Float64,Int64}(; skip_points, buffer_size, save_file,
                                        auto_notify=false)
    return ActiveHyperedgeCount(log, size)
end

mutable struct MotifCount <: AbstractStepMeasurement
    log::MeasurementLog{Float64,Int64}
    label::Label
end

function MotifCount(label::Label; skip_points::Int64=1, buffer_size::Int64=0,
                    save_folder::Union{Nothing,String})
    save_file = _create_save_file(save_folder, "motif_count", "$label")
    log = MeasurementLog{Float64,Int64}(; skip_points, buffer_size, save_file,
                                        auto_notify=false)
    return MotifCount(log, label)
end

mutable struct FakeDiffEq <: AbstractStepMeasurement
    log::MeasurementLog{Float64,Float64}
    label::Label
end

function FakeDiffEq(label::Label; skip_points::Int64=1, buffer_size::Int64=0,
                    save_folder::Union{Nothing,String})
    save_file = _create_save_file(save_folder, "fake_diff_eq", "$label")
    log = MeasurementLog{Float64,Float64}(; skip_points, buffer_size, save_file,
                                          auto_notify=false)
    return FakeDiffEq(log, label)
end

# ====================================================================================
# ------------------------------- Run Measurements -----------------------------------

"""
Measures the time that the system needs to deplete all active hyperedges. 
"""
struct ActiveLifetime <: AbstractRunMeasurement
    log::MeasurementLog{Int64,Float64}
end

function ActiveLifetime(; save_folder::Union{Nothing,String})
    save_file = _create_save_file(save_folder, "active_lifetime")
    log = MeasurementLog{Int64,Float64}(; save_file=save_file)
    return ActiveLifetime(log)
end

"""
Measures the magnetization in the end of the simulation. 
"""
struct FinalMagnetization <: AbstractRunMeasurement
    log::MeasurementLog{Int64,Int64}
end

function FinalMagnetization(; save_folder::Union{Nothing,String})
    save_file = _create_save_file(save_folder, "final_magnetization")
    log = MeasurementLog{Int64,Float64}(; save_file=save_file)
    return FinalMagnetization(log)
end

struct AvgHyperedgeCount <: AbstractRunMeasurement
    log::MeasurementLog{Int64,Float64}
    label::Int64
end

function AvgHyperedgeCount(size::Int64; save_folder::Union{Nothing,String})
    save_file = _create_save_file(save_folder, "avg_hyperedge_count")
    return AvgHyperedgeCount(MeasurementLog{Int64,Float64}(; save_file=save_file), size)
end

struct SlowManifoldFit <: AbstractRunMeasurement
    log::MeasurementLog{Int64,NTuple{3,Float64}}
    label::Int64
end

function SlowManifoldFit(size::Int64; save_folder::Union{Nothing,String})
    save_file = _create_save_file(save_folder, "slow_manifold_fit")
    log = MeasurementLog{Int64,NTuple{3,Float64}}(; save_file=save_file)
    return SlowManifoldFit(log, size)
end

"""
_snake_case(str:S) where S <: AbstractString

Helper function to convert type names in CamelCase to property names in snake_case.
    
    Copied from: https://stackoverflow.com/questions/70007955/julia-implementation-for-converting-string-to-snake-case-camelcase
    """
function _snake_case(str::S) where {S<:AbstractString}
    wordpat = r"
    ^[a-z]+ |                  #match initial lower case part
    [A-Z][a-z]+ |              #match Words Like This
    \d*([A-Z](?=[A-Z]|$))+ |   #match ABBREV 30MW 
    \d+                        #match 1234 (numbers without units)
    "x

    smartlower(word) = any(islowercase, word) ? lowercase(word) : word
    words = [smartlower(m.match) for m in eachmatch(wordpat, str)]

    return join(words, "_")
end

function set_save_file!(meas::AbstractMeasurement, save_folder::Union{Nothing,String})
    # don't do anything if the folder didn't change
    old_path = meas.log.save_file
    if !isnothing(old_path) && save_folder == splitdir(old_path)[1]
        return meas
    end

    meas_str = _snake_case("$(typeof(meas))")
    label = :label in fieldnames(typeof(meas)) ? getfield(meas, :label) : ""
    save_file = _create_save_file(save_folder, meas_str, label)
    meas.log.save_file = save_file
    return meas
end

function _create_save_file(save_folder::Union{Nothing,String}, meas::String,
                           label="")
    if isnothing(save_folder)
        return nothing
    end
    filename = if (label != "")
        "$(meas)_$label.csv"
    else
        "$meas.csv"
    end
    save_file = joinpath(save_folder, filename)
    mkpath(save_folder)
    open(save_file, "w") do io
        return write(io, "id, index, value\n")
    end
    return save_file
end