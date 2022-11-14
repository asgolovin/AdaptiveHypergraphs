using Parameters
using DrWatson
using JSON3

export InputParams, NetworkParams, ModelParams, VisualizationParams, BatchParams, save_json

@with_kw struct NetworkParams
    num_nodes::Union{Int64,Vector{Int64}} = 100
    infected_prob::Union{Float64,Vector{Float64}} = 0.5
    num_hyperedges::Union{Tuple,Vector{<:Tuple}} = (100, 10)
end

@with_kw struct ModelParams
    is_discrete::Union{Bool,Vector{Bool}} = true
    adaptivity_rule::Union{AdaptivityRule,Vector{<:AdaptivityRule}} = RewireToRandom()
    propagation_rule::Union{PropagationRule,Vector{<:PropagationRule}} = MajorityVoting()
    num_time_steps::Union{Int64,Vector{Int64}} = 500
    adaptivity_prob::Union{Float64,Vector{Float64}} = 0.5
    propagation_rate::Union{Float64,Vector{Float64}} = 1.0
    adaptivity_rate::Union{Float64,Vector{Float64}} = 1.0
end

@with_kw struct VisualizationParams
    skip_points::Int64 = 1
    buffer_size::Int64 = 100
    node_colormap::Symbol = :RdYlGn_6
    hyperedge_colormap::Symbol = :thermal
end

@with_kw struct BatchParams
    record_video::Bool = false
    batch_size::Int64 = 10
    # turns on a prompt if the data should be saved. The prompt is a bit annoying, 
    # so it can be turned off completely with this option. 
    prompt_for_save::Bool = false
end

struct InputParams
    network_params::NetworkParams
    model_params::ModelParams
    visualization_params::VisualizationParams
    batch_params::BatchParams
end

# Find all parameters which are given as vectors
function _find_sweeps(params::InputParams)
    param_groups = fieldnames(InputParams)
    sweeps = []
    for group in [params.network_params, params.model_params]
        for param in fieldnames(typeof(group))
            value = getproperty(group, param)
            if typeof(value) <: Vector
                push!(sweeps, param)
            end
        end
    end
    return sweeps
end

"""
    flatten(params::InputParams)

Return a vector with all combinations of all sweeped params. 
"""
function flatten(params::InputParams)
    sweeps = _find_sweeps(params)
    @show sweeps
    if length(sweeps) == 0
        return [params]
    end

    param_vector = InputParams[]

    for sweep in sweeps
        for value in getproperty(params.model_params, sweep)
            param_copy = params
            setproperty!(param_copy.model_params, sweep, value)
            param_vector.append(param_copy)
        end
    end
    return param_vector
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