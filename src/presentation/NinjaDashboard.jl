using DrWatson

export NinjaDashboard, PANEL_DEPENDENCIES

#! format: off
"""
A list of dependencies of Panels on Measurements. 
"""
PANEL_DEPENDENCIES = Dict{Symbol, Vector{DataType}}(
        :HypergraphPanel               => [],
        :StateDistPanel                => [StateCount],
        :HyperedgeDistPanel            => [HyperedgeCount, AvgHyperedgeCount],
        :ActiveHyperedgeDistPanel      => [ActiveHyperedgeCount],
        :FirstOrderMotifCountPanel     => [MotifCount],
        :MomentClosurePanel            => [MotifCount, StateCount],
        :SecondOrderMotifCountPanel    => [MotifCount],
        :FakeDiffEqPanel               => [FakeDiffEq, MotifCount, StateCount],
        :ActiveRatioPanel              => [ActiveHyperedgeCount],
        :SlowManifoldPanel             => [StateCount, ActiveHyperedgeCount, SlowManifoldFit],
        :ActiveLifetimePanel           => [ActiveLifetime],
        :FinalMagnetizationPanel       => [FinalMagnetization],
        :AvgHyperedgeCountPanel        => [AvgHyperedgeCount, SlowManifoldFit])
#! format: on

abstract type AbstractDashboard end

"""
    NinjaDashboard <: AbstractDashboard

A dashboard, but it's so stealthy that you will never see it. 
It does not use any graphical libraries. 
It never displays anything at all, just computes the data in the darkness of the night. 
"""
struct NinjaDashboard <: AbstractDashboard
    mo::ModelObservable
    measurement_types::Vector{DataType}
end

"""
    NinjaDashboard(model::AbstractModel)

# Arguments
- `model::AbstractModel` - the model on which the dashboard is based on. 
"""
function NinjaDashboard(model::AbstractModel, vparams::VisualizationParams;
                        save_folder::Union{Nothing,String})
    # collect a list of the measurements on which the panels depend on
    measurements = Vector{DataType}()
    for panel in vparams.panels
        append!(measurements, PANEL_DEPENDENCIES[panel])
    end
    measurement_types = unique(measurements)

    # throw out the run measurements because we don't have all the data
    measurement_types = filter!(x -> x <: AbstractStepMeasurement, measurement_types)

    mo = ModelObservable(model, measurement_types;
                         skip_points=vparams.skip_points,
                         buffer_size=vparams.buffer_size,
                         write_to_observables=false,
                         save_folder=save_folder)
    return NinjaDashboard(mo, measurement_types)
end

"""
    run!(dashboard::AbstractDashboard, num_steps::Int64)

Run the simulation for the amount of time set by `duration` or until the hypergraph runs out of active hyperedges.
"""
function run!(dashboard::AbstractDashboard, duration::Float64)
    mo = dashboard.mo
    i = mo.num_steps
    active_lifetime = duration

    while mo.time <= duration
        step!(mo)
        i += 1
        num_active_hyperedges = get_num_active_hyperedges(mo.network[])

        # a hack to say that the type has to be Dashboard without mentioning the Dashboard
        if !(typeof(dashboard) <: NinjaDashboard)
            if i % mo.buffer_size == 0 || num_active_hyperedges == 0
                for panel in dashboard.panels
                    set_lims!(panel)
                end
            end
        end

        # stop the simulation early if we run out of active hyperdeges
        if num_active_hyperedges == 0
            active_lifetime = i
            break
        end
    end
    record_measurements!(mo, :step)
    flush_buffers!(mo)
    notify(mo)

    return dashboard
end

"""
    reset!(dashboard::NinjaDashboard, model::AbstractModel)

Reset the dashboard to run the next simulation from the batch. 
The observables in the ModelObservable are reset to track the new data from `model`.
"""
function reset!(dashboard::NinjaDashboard, model::AbstractModel,
                save_folder::Union{Nothing,String})
    # reset observables
    rebind_model!(dashboard.mo, model, save_folder)
    return dashboard
end