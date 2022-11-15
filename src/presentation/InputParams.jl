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

@with_kw struct InputParams
    network_params::NetworkParams = NetworkParams()
    model_params::ModelParams = ModelParams()
    visualization_params::VisualizationParams = VisualizationParams()
    batch_params::BatchParams = BatchParams()
end

"""
    expand(params::InputParams)

Expand the params into a vector of params which contains all possible combinations 
over all ranged params.

For example, say, we want to run several simulations over the following ranges of parameters:

```
    params.example_param = [1, 2, 3]
    params.second_param = ["a", "b"]
```

The function produces a list of params with a Cartesian product of both sets:
{{1, "a"}, {2, "a"}, {3, "a"}, {1, "b"}, {2, "b"}, {3, "b"}}. Non-vectored params are left unchanged. 

For a param to be expanded, it has to have type Vector; other iterables are not supported. 
Furthermore, only vectors in NetworkParams and ModelParams are supported, since it doesn't make much sense to change visualization and batch params between batches.
"""
function expand(params::InputParams)
    vparams, bparams = params.visualization_params, params.batch_params

    # convert the InputParams struct into a dict
    nparams_dict = Dict(key => getfield(params.network_params, key)
                        for key in fieldnames(NetworkParams))
    mparams_dict = Dict(key => getfield(params.model_params, key)
                        for key in fieldnames(ModelParams))
    params_dict = merge(nparams_dict, mparams_dict)

    # expand the dict into a vector of dicts
    flattened_dicts = dict_list(params_dict)

    # convert the dict back to InputParams
    flattened_params = []
    for d in flattened_dicts
        # filter out only the relevant params for each param type
        nargs = filter(x -> x.first in fieldnames(NetworkParams), d)
        margs = filter(x -> x.first in fieldnames(ModelParams), d)
        nparams = NetworkParams(; nargs...)
        nparams.num_nodes
        mparams = ModelParams(; margs...)
        push!(flattened_params, InputParams(nparams, mparams, vparams, bparams))
    end
    return flattened_params
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