using Statistics, Parameters, Polynomials

export ModelObservable, step!, flush_buffers!, notify, record_measurements!,
       rebind_model!, clear!

#! format: off
"""
A list of dependencies of Measurements on other Measurements. 
"""
MEASUREMENT_DEPENDENCIES = Dict{DataType, Vector{DataType}}(
        StateCount              => [],
        HyperedgeCount          => [],
        ActiveHyperedgeCount    => [],
        MotifCount              => [],
        FakeDiffEq              => [],
        ActiveLifetime          => [],
        AvgHyperedgeCount       => [HyperedgeCount, ActiveHyperedgeCount],
        FinalMagnetization      => [],
        SlowManifoldFit         => [StateCount, ActiveHyperedgeCount])
#! format: on

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
@with_kw mutable struct ModelObservable
    time::Float64
    num_steps::Int64
    skip_points::Int64
    buffer_size::Int64
    model::Observable{AbstractModel}
    network::Observable{HyperNetwork}
    state_count::Vector{StateCount} = StateCount[]
    hyperedge_count::Vector{HyperedgeCount} = HyperedgeCount[]
    active_hyperedge_count::Vector{ActiveHyperedgeCount} = ActiveHyperedgeCount[]
    motif_count::Vector{MotifCount} = MotifCount[]
    fake_diff_eq::Vector{FakeDiffEq} = FakeDiffEq[]
    active_lifetime::Vector{ActiveLifetime} = ActiveLifetime[]
    final_magnetization::Vector{FinalMagnetization} = FinalMagnetization[]
    avg_hyperedge_count::Vector{AvgHyperedgeCount} = AvgHyperedgeCount[]
    slow_manifold_fit::Vector{SlowManifoldFit} = SlowManifoldFit[]
end

function ModelObservable(model::AbstractModel, measurement_types::Vector{DataType};
                         skip_points=1, buffer_size=1,
                         save_folder::Union{Nothing,String}=nothing)
    time = 0.0
    num_steps = 0

    max_size = get_max_size(model.network)

    # the required measurements from measurement_types might depend on other measurements. 
    # We need to add them to the list.
    buffer = copy(measurement_types)
    measurement_types_expanded = []
    while length(buffer) > 0
        mtype = pop!(buffer)
        if !(mtype in measurement_types_expanded)
            push!(measurement_types_expanded, mtype)
        end
        for dependency in MEASUREMENT_DEPENDENCIES[mtype]
            push!(buffer, dependency)
        end
    end

    # instantiate only the required measurements
    measurements = Dict()
    log_params = Dict(:skip_points => skip_points,
                      :buffer_size => buffer_size)
    for type in measurement_types_expanded
        type_str = _snake_case("$type")
        sym = Symbol(type_str)
        if typeof(save_folder) <: String
            save_folder = joinpath([save_folder, type_str])
        end
        if type <: StateCount
            measurements[sym] = [StateCount(state; save_folder=save_folder, log_params...)
                                 for state in instances(State)]
        elseif type <: HyperedgeCount || type <: ActiveHyperedgeCount
            measurements[sym] = [type(size; save_folder=save_folder, log_params...)
                                 for size in 2:max_size]
        elseif type <: MotifCount
            measurements[sym] = [type(label; save_folder=save_folder, log_params...)
                                 for label in all_labels(max_size)
                                 if order(label) > 0]
        elseif type <: FakeDiffEq
            measurements[sym] = [type(label; save_folder=save_folder, log_params...)
                                 for label in all_labels(max_size)
                                 if order(label) == 1 && label.left[A] > 0 &&
                                    label.left[B] > 0]
        elseif type <: AvgHyperedgeCount || type <: SlowManifoldFit
            measurements[sym] = [type(size; save_folder=save_folder) for size in 2:max_size]
        else
            measurements[sym] = [type(; save_folder=save_folder)]
        end
    end

    return ModelObservable(; time=time,
                           num_steps=num_steps,
                           skip_points=skip_points,
                           buffer_size=buffer_size,
                           model=Observable(model),
                           network=Observable(model.network),
                           measurements...)
end

function Base.show(io::IO, mo::ModelObservable)
    return print("ModelObservable(model::$(typeof(mo.model[])))")
end

function Base.show(io::IO, ::MIME"text/plain", mo::ModelObservable)
    println(io, "ModelObservable of a $(typeof(mo.model[]))")
    println(io, "  time: $(mo.time)")
    println(io, "  num_steps: $(mo.num_steps)")
    step_measurements = join(mo.step_measurements, ", ")
    run_measurements = join(mo.run_measurements, ", ")
    batch_measurements = join(mo.batch_measurements, ", ")
    println(io, "  step measurements: [$step_measurements]")
    println(io, "  run measurements: [$run_measurements]")
    return println(io, "  batch measurements: [$batch_measurements]")
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
    elseif sym === :batch_measurements
        return filter(x -> typeof(x) <: AbstractBatchMeasurement, obj.measurements)
    else
        return getfield(obj, sym)
    end
end

"""
    step!(mo::ModelObservable)
    
Progress the model one time step forward and update the history. 
"""
function step!(mo::ModelObservable)
    network_changed, Δt = step!(mo.model[])
    mo.time += Δt
    mo.num_steps += 1
    network_changed && notify(mo.model)
    record_measurements!(mo, :step)
    if mo.num_steps % mo.buffer_size == 0
        sleep(0.001)
        notify(mo)
        println(mo.time)
    end
    return mo
end

"""
    write(mo::ModelObservable)

Save the new data collected since the last save to csv files.
"""
function checkpoint_save(mo::ModelObservable)
    # TODO
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

context can be either `:step`, `:run` or `:batch`.
"""
function record_measurements!(mo::ModelObservable, context::Symbol)
    if context == :step
        measurements = mo.step_measurements
    elseif context == :run
        measurements = mo.run_measurements
    elseif context == :batch
        measurements = mo.batch_measurements
    else
        raise(ArgumentError("context should be either :step, :run or :batch."))
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
function rebind_model!(mo::ModelObservable, model::AbstractModel,
                       save_folder::Union{Nothing,String})
    clear!(mo)
    mo.time = 0.0
    mo.num_steps = 0
    mo.model[] = model
    mo.network[] = model.network

    for measurement in mo.measurements
        set_save_file!(measurement, save_folder)
    end

    return mo
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

function record_measurement!(mo::ModelObservable, measurement::MotifCount)
    label = measurement.label
    count = get_motif_count(mo.network[])[label]
    record!(measurement.log, mo.time, count)
    return measurement
end

function record_measurement!(mo::ModelObservable, measurement::FakeDiffEq)
    label = measurement.label
    if mo.num_steps == 1
        value = Float64(get_motif_count(mo.network[])[label])
        record!(measurement.log, mo.time, value)
        return measurement
    end

    k = label.left_total[A]
    h = label.left_total[B]
    p = mo.model[].adaptivity_prob
    max_size = get_max_size(mo.network[])

    hyperedge_dist = get_hyperedge_dist(mo.network[])
    # number of hyperedges which can be rewired to, i.e., those which are not of max size
    num_small_hyperedges = sum([x.second for x in hyperedge_dist if x.first != max_size])
    num_nodes = get_num_nodes(mo.network[])
    # total number of candidates which can be rewired to
    num_candidates = num_small_hyperedges + num_nodes

    motif_dict = get_motif_count(mo.network[])

    # the difference between the value at t and t+Δt
    update = 0

    # PROPAGATION

    prop_update = 0.0

    for n in 1:max_size, m in 1:(max_size - n)
        if n + m == 1
            continue
        end
        prop_update += n / (n + m) * (motif_dict[Label("[A$(n) B$(m-1)|B|A$(k-1)B$(h)]")] -
                                      motif_dict[Label("[A$(n) B$(m-1)|B|A$(k)B$(h-1)]")])
        prop_update += m / (n + m) * (motif_dict[Label("[A$(n-1) B$(m)|A|A$(k)B$(h-1)]")] -
                                      motif_dict[Label("[A$(n-1) B$(m)|A|A$(k-1)B$(h)]")])
    end
    # symmetric terms
    #! format: off
    prop_update += (k - 1) / (k + h) * motif_dict[Label("[A$(k-1)B$(h)|B|A$(k-1)B$(h)]")] -
                         k / (k + h) * motif_dict[Label("[A$(k)B$(h-1)|B|A$(k)B$(h-1)]")] +
                   (h - 1) / (k + h) * motif_dict[Label("[A$(k)B$(h-1)|A|A$(k)B$(h-1)]")] -
                         h / (k + h) * motif_dict[Label("[A$(k-1)B$(h)|A|A$(k-1)B$(h)]")]

    
    # ADAPTIVITY
    adapt_update = 0.0

    for n in 1:max_size, m in 1:(max_size - n)
        if n + m == 1
            continue
        end
        AnBm = motif_dict[Label("[A$(n)B$(m)]")]

        adapt_update += n / (n + m) * motif_dict[Label("[A$(k-1)B$(h)]")] * AnBm / num_candidates +
                        m / (n + m) * motif_dict[Label("[A$(k)B$(h-1)]")] * AnBm / num_candidates -
                        motif_dict[Label("[A$(k-1)B$(h)]")] * AnBm / num_candidates
    end

    if k + h < max_size
        adapt_update += (k + 1) / (k + h + 1) * motif_dict[Label("[A$(k+1)B$(h)]")]
        adapt_update += (h + 1) / (k + h + 1) * motif_dict[Label("[A$(k)B$(h+1)]")]
    end

    #! format: on
    last_meas = measurement.log.last_value
    num_hyperedges = get_num_hyperedges(mo.network[])
    update = ((1 - p) * prop_update + p * adapt_update - last_meas) / num_hyperedges

    value = last_meas + update

    record!(measurement.log, mo.time, value)

    return measurement
end

function record_measurement!(mo::ModelObservable, measurement::ActiveLifetime)
    record!(measurement.log, mo.time)
    return measurement
end

function record_measurement!(mo::ModelObservable, measurement::FinalMagnetization)
    state_count = get_state_count(mo.network[])
    magnetization = state_count[A] - state_count[B]
    record!(measurement.log, magnetization)
    return measurement
end

function record_measurement!(mo::ModelObservable, measurement::AvgHyperedgeCount)
    size = measurement.label
    hyperedge_timeseries = mo.hyperedge_count[size - 1].values[]
    skip = length(hyperedge_timeseries) ÷ 10
    avg_hyperedge_count = mean(hyperedge_timeseries[skip:end])
    record!(measurement.log, avg_hyperedge_count)
    return measurement
end

function record_measurement!(mo::ModelObservable, measurement::SlowManifoldFit)
    size = measurement.label
    x = mo.state_count[1].values[]
    y = mo.active_hyperedge_count[size - 1].values[]
    skip = length(x) ÷ 10
    x = x[skip:end]
    y = y[skip:end]

    f = Polynomials.fit(x, y, 2) # polynomial fit of degree 2
    a, b, c = coeffs(f) # f(x) = a + bx + cx^2
    x_peak = -b / (2 * c)
    peak = a + b * x_peak + c * x_peak^2
    println("size: $size, peak at $peak")
    record!(measurement.log, (a, b, c))
    return measurement
end