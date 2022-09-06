using Parameters
using DrWatson
using JSON3

export InputParams, NetworkParams, ModelParams, VisualizationParams, BatchParams, save_json

@with_kw struct NetworkParams
    num_nodes::Integer = 100
    infected_prob::Real = 0.5
    num_hyperedges::Tuple = (100, 10)
end

@with_kw struct ModelParams
    is_discrete::Bool = true
    adaptivity_rule::AdaptivityRule = RewireToRandom()
    propagation_rule::PropagationRule = MajorityVoting()
    num_time_steps::Integer = 500
    propagation_prob::Real = 0.5
    propagation_rate::Real = 1.0
    adaptivity_rate::Real = 1.0
end

@with_kw struct VisualizationParams
    dashboard_params::NamedTuple = ()
    steps_per_update::Integer = 10
end

@with_kw struct BatchParams
    record_video::Bool = false
    batch_size::Integer = 10
end

struct InputParams
    network_params::NetworkParams
    model_params::ModelParams
    visualization_params::VisualizationParams
    batch_params::BatchParams
end

function save_json(io::IO, params::InputParams)
    return JSON3.pretty(io, save_json(params))
end

function save_json(params::InputParams)
    json_dict = Dict{Symbol,Any}()
    json_dict[:network_params] = struct2dict(params.network_params)
    json_dict[:model_params] = struct2dict(params.model_params)
    json_dict[:model_params][:propagation_rule] = nameof(typeof(json_dict[:model_params][:propagation_rule]))
    json_dict[:model_params][:adaptivity_rule] = nameof(typeof(json_dict[:model_params][:adaptivity_rule]))
    json_dict[:visualization_params] = struct2dict(params.visualization_params)
    json_dict[:batch_params] = struct2dict(params.batch_params)

    return JSON3.write(json_dict)
end