using Statistics, Parameters

export ModelObservable, step!, flush_buffers!, notify, record_measurements!,
       rebind_model!, clear!

"""
    ModelObservable{M <: AbstarctModel}

Contains all observables (both as in physical observables and as in Makie Observable class) 
needed to plot the evolution of the model. 

This class be methaphorically thought of as a scientist who runs some experiments 
(i.e., pushes a button to evolve the model one step forward), observes the results 
and puts a new marker into a plot. ModelObservable provides its own interface to the 
step! function which evolves the model and records the results. 

ModelObservable gets a list of required measurements from Dashboard and collects data
on those measurements only. 

To add a new measurement:
    1. create a new struct in Measuremnts.jl either of type AbstractStepMeasurement or AbstractRunMeasurement
    2. add a `record_measurement!(mo::ModelObservable, measurement::YourNewMeasurement)` function.
    3. make the measurement a property of ModelObservable. The name of the field should be the name of the struct converted to snake_case as returned by `_snake_case(YourNewMeasurement)`
    4. it might be necessary to add a case distinction to create the measurement in the constructor of `ModelObservable`. 
"""
@with_kw mutable struct ModelObservable{M<:AbstractModel}
    time::Int64
    buffer_size::Int64
    model::Observable{M}
    network::Observable{HyperNetwork}
    state_count::Vector{StateCount} = StateCount[]
    hyperedge_count::Vector{HyperedgeCount} = HyperedgeCount[]
    active_hyperedge_count::Vector{ActiveHyperedgeCount} = ActiveHyperedgeCount[]
    active_lifetime::Vector{ActiveLifetime} = ActiveLifetime[]
    final_magnetization::Vector{FinalMagnetization} = FinalMagnetization[]
    avg_hyperedge_count::Vector{AvgHyperedgeCount} = AvgHyperedgeCount[]
end

function ModelObservable{M}(model::M, measurement_types::Vector{DataType};
                            skip_points=1,
                            buffer_size=1) where {M<:AbstractModel}
    time = 0

    max_size = get_max_hyperedge_size(model.network)
    # instantiate only the required measurements
    measurements = Dict()
    log_params = Dict(:skip_points => skip_points,
                      :buffer_size => buffer_size)
    for type in measurement_types
        sym = Symbol(_snake_case("$type"))
        if type <: StateCount
            measurements[sym] = [StateCount(state; log_params...)
                                 for state in instances(State)]
        elseif type <: HyperedgeCount || type <: ActiveHyperedgeCount
            measurements[sym] = [type(size; log_params...)
                                 for size in 2:max_size]
        elseif type <: AvgHyperedgeCount
            measurements[sym] = [type(size) for size in 2:max_size]
        else
            measurements[sym] = [type()]
        end
    end

    return ModelObservable(; time=time,
                           buffer_size=buffer_size,
                           model=Observable(model),
                           network=Observable(model.network),
                           measurements...)
end

function Base.getproperty(obj::ModelObservable, sym::Symbol)
    if sym === :measurements
        names = fieldnames(ModelObservable)
        values = [getproperty(obj, name) for name in names]
        filter!(x -> typeof(x) <: Vector{<:AbstractMeasurement}, values)
        return vcat(values...)
    elseif sym === :step_measurements
        return filter(x -> typeof(x) <: AbstractStepMeasurement, obj.measurements)
    elseif sym === :run_measurements
        return filter(x -> typeof(x) <: AbstractRunMeasurement, obj.measurements)
    else
        return getfield(obj, sym)
    end
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

"""
    step!(mo::ModelObservable)
    
Progress the model one time step forward and update the history. 
"""
function step!(mo::ModelObservable)
    network_changed = step!(mo.model[])
    network_changed && notify(mo.model)
    record_measurements!(mo, :step)
    if mo.time % mo.buffer_size == 0
        sleep(0.001)
        notify(mo)
    end
    mo.time += 1
    return mo
end

"""
    flush_buffers!(mo::ModelObservable)

Write all values in buffers of logs to the observables. 

Used to flush the buffers "by hand" in the end of the simulation. 
"""
function flush_buffers!(mo::ModelObservable)
    for measurement in mo.step_measurements
        flush_buffers!(measurement.log)
    end
    for measurement in mo.step_measurements
        notify(measurement.log)
    end
    return mo
end

"""
    notify(mo::ModelObservable)

Notify all observables in the ModelObservable.
"""
function GLMakie.notify(mo::ModelObservable)
    notify(mo.network)
    for measurement in mo.measurements
        notify(measurement.log)
    end
    return nothing
end

"""
    record_measurements!(mo::ModelObservable)

Record measurements which fit the corresponding context. 

context can be either `:step` or `:run`.
"""
function record_measurements!(mo::ModelObservable, context::Symbol)
    if context == :step
        measurements = mo.step_measurements
    elseif context == :run
        measurements = mo.run_measurements
    else
        raise(ArgumentError("context should be either :step or :run."))
    end

    for measurement in measurements
        record_measurement!(mo, measurement)
    end
    return mo
end

"""
    rebind_model!(mo::ModelObservable, model::AbstractModel)

Rebind the observables to track the new `model`.
"""
function rebind_model!(mo::ModelObservable, model::AbstractModel)
    clear!(mo)
    mo.time = 0
    mo.model[] = model
    mo.network[] = model.network
    return record_measurements!(mo, :step)
end

"""
    clear!(mo::ModelObservable)

Clear the history buffers in the ModelObservable.
"""
function clear!(mo::ModelObservable)
    for measurement in mo.step_measurements
        clear!(measurement.log)
    end
    return mo
end

function record_measurement!(mo::ModelObservable, measurement::StateCount)
    state_count = get_state_count(mo.network[])
    state = measurement.label
    record!(measurement.log, mo.time, state_count[state])
    return measurement
end

function record_measurement!(mo::ModelObservable, measurement::HyperedgeCount)
    hyperedge_count = get_hyperedge_dist(mo.network[])
    size = measurement.label
    record!(measurement.log, mo.time, hyperedge_count[size])
    return measurement
end

function record_measurement!(mo::ModelObservable, measurement::ActiveHyperedgeCount)
    size = measurement.label
    count = get_num_active_hyperedges(mo.network[], size)
    record!(measurement.log, mo.time, count)
    return measurement
end

function record_measurement!(mo::ModelObservable, measurement::ActiveLifetime)
    record!(measurement.log, mo.time)
    return measurement
end

function record_measurement!(mo::ModelObservable, measurement::FinalMagnetization)
    state_count = get_state_count(mo.network[])
    magnetization = state_count[S] - state_count[I]
    record!(measurement.log, magnetization)
    return measurement
end

function record_measurement!(mo::ModelObservable, measurement::AvgHyperedgeCount)
    size = measurement.label
    avg_hyperedge_count = mean(mo.hyperedge_count[size - 1].log.values[])
    avg_active_count = mean(mo.active_hyperedge_count[size - 1].log.values[])
    record!(measurement.log, (total=avg_hyperedge_count, active=avg_active_count))
    return measurement
end