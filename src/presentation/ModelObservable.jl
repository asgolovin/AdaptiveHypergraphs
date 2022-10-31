using Statistics

export ModelObservable, step!, flush_buffers!, record_time_series, record_active_lifetime!,
       rebind_model!, clear!

"""
    ModelObservable{M <: AbstarctModel}

Contains all observables (both as in physical observables and as in Makie Observable class) 
needed to plot the evolution of the model. 

This class be methaphorically thought of as a scientist who runs some experiments 
(i.e., pushes a button to evolve the model one step forward), observes the results 
and puts a new marker into a plot. ModelObservable provides its own interface to the 
step! function which evolves the model and records the results. 
"""
mutable struct ModelObservable{M<:AbstractModel}
    time::Int64
    model::Observable{M}
    network::Observable{HyperNetwork}
    measurements::Vector

    function ModelObservable{M}(model::M, skip_points=1,
                                buffer_size=1) where {M<:AbstractModel}
        time = 0
        measurements = []
        state_series = StateCount(skip_points, buffer_size)
        push!(measurements, state_series)

        max_size = get_max_hyperedge_size(model.network)
        hyperedge_series = HyperedgeCount(max_size, skip_points, buffer_size)
        push!(measurements, hyperedge_series)

        active_hyperedges_series = ActiveHyperedgeCount(max_size, skip_points, buffer_size)
        push!(measurements, active_hyperedges_series)

        final_hyperedge_dist = FinalHyperedgeDist()
        push!(measurements, final_hyperedge_dist)

        active_lifetime = ActiveLifetime()
        push!(measurements, active_lifetime)

        final_magnetization = FinalMagnetization()
        push!(measurements, final_magnetization)

        mo = new{M}(time,
                    Observable(model),
                    Observable(model.network),
                    measurements)
        return record_measurements!(mo, :step)
    end
end

"""
    step!(mo::ModelObservable)
    
Progress the model one time step forward and update the history. 
"""
function step!(mo::ModelObservable)
    network_changed = step!(mo.model[])
    network_changed && notify(mo.model)
    record_measurements!(mo, :step)
    mo.time += 1
    return mo
end

"""
record_measurements!(mo::ModelObservable)

Record measurements which fit the corresponding context. 

context can be either :step or :run
"""
function record_measurements!(mo::ModelObservable, context::Symbol)
    if context == :step
        measurement_type = AbstractStepMeasurement
    elseif context == :run
        measurement_type = AbstractRunMeasurement
    else
        raise(ArgumentError("context should be either :step or :run."))
    end

    for measurement in mo.measurements
        if typeof(measurement) <: measurement_type
            record_measurement!(mo, measurement)
        end
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
    for series in _get_all_series(mo)
        clear!(series)
    end
    return mo
end

function record_measurement!(mo::ModelObservable, measurement::StateCount)
    state_count = get_state_count(mo.network[])
    for state in instances(State)
        record!(measurement.log[state], mo.time, state_count[state])
    end
    return measurement
end

function record_measurement!(mo::ModelObservable, measurement::HyperedgeCount)
    hyperedge_count = get_hyperedge_dist(mo.network[])
    for size in 2:get_max_hyperedge_size(mo.network[])
        record!(measurement.log[size], mo.time, hyperedge_count[size])
    end
    return measurement
end

function record_measurement!(mo::ModelObservable, measurement::ActiveHyperedgeCount)
    for size in 2:get_max_hyperedge_size(mo.network[])
        count = get_num_active_hyperedges(mo.network[], size)
        record!(measurement.log[size], mo.time, count)
    end
    return measurement
end

function record_measurement!(mo::ModelObservable, measurement::ActiveLifetime)
    record!(measurement.log, mo.time)
    return measurement
end

function record_measurement!(mo::ModelObservable, measurement::FinalMagnetization)
    state_count = get_state_count(mo.network[])
    magnetization = state_count[S] - state_count[I]
    record!(measurement, magnetization)
    return measurement
end

function record_measurement!(mo::ModelObservable, measurement::FinalHyperedgeDist)
    # We want to calculate the average value of the hyperdeges after the system has stabilized. 
    # Idea: compute the std of the time series starting from t to the end of the simulation 
    std_series = Dict{Int64,Vector{Real}}()
    for hyperedge_count in mo.hyperedge_series
        size = hyperedge_count.size
        std_series[size] = []
        for t in 1:(length(hyperedge_count.values[]) - 1)
            push!(std_series[size], Statistics.std(hyperedge_count.values[][t:end]))
        end
    end
    record!(measurement, collect(values(std_series)))
    return measurement
end