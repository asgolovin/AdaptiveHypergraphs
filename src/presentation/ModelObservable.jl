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
struct ModelObservable{M<:AbstractModel}
    model::Observable{M}
    network::Observable{HyperNetwork}
    state_series::Vector{StateCount}
    hyperedge_series::Vector{HyperedgeCount}
    active_hyperedges_series::ActiveHyperedgeCount
    active_lifetime::ActiveLifetime
    final_magnetization::FinalMagnetization

    function ModelObservable{M}(model::M) where {M<:AbstractModel}
        state_series = [StateCount(model.network, state) for state in instances(State)]
        max_size = get_max_hyperedge_size(model.network)
        hyperedge_series = [HyperedgeCount(model.network, size) for size in 2:max_size]
        active_hyperedges_series = ActiveHyperedgeCount(model.network)
        active_lifetime = ActiveLifetime()
        final_magnetization = FinalMagnetization()

        mo = new{M}(Observable(model),
                    Observable(model.network),
                    state_series,
                    hyperedge_series,
                    active_hyperedges_series,
                    active_lifetime,
                    final_magnetization)
        return record_time_series!(mo)
    end
end

function _get_all_series(mo::ModelObservable)
    return [mo.state_series; mo.hyperedge_series; mo.active_hyperedges_series]
end

"""
    step!(mo::ModelObservable)
    
Progress the model one time step forward and update the history. 
"""
function step!(mo::ModelObservable)
    network_changed = step!(mo.model[])
    network_changed && notify(mo.model)
    return record_time_series!(mo)
end

function flush_buffers!(mo::ModelObservable)
    for series in _get_all_series(mo)
        flush_buffers!(series)
    end
    for series in _get_all_series(mo)
        notify(series.observable)
    end
    return mo
end

"""
    record_time_series!(mo::ModelObservable)

Push the current distribution of states from the model into the history buffer vectors. 
"""
function record_time_series!(mo::ModelObservable)
    for series in _get_all_series(mo)
        record_measurement!(series)
    end
    return mo
end

function record_active_lifetime!(mo::ModelObservable, time)
    record_measurement!(mo.active_lifetime, time)
    return mo
end

"""
    rebind_model!(mo::ModelObservable, model::AbstractModel)

Rebind the observables to track the new `model`.
"""
function rebind_model!(mo::ModelObservable, model::AbstractModel)
    clear!(mo)
    mo.model[] = model
    mo.network[] = model.network
    for series in _get_all_series(mo)
        series.network = model.network
    end
    return record_time_series!(mo)
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