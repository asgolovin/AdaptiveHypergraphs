using Parameters
using DrWatson
using JSON3

export InputParams, NetworkParams, ModelParams, VisualizationParams, BatchParams, save_json,
       load_params

@with_kw mutable struct NetworkParams
    # Number of nodes in the hypergraph
    num_nodes::Union{Int64,Vector{Int64}} = 100
    # Probability that a given node is in state A
    state_A_prob::Union{Float64,Vector{Float64}} = 0.5
    # Number of hyperedges of every size in the hypergraph. For example, (100, 10) means that there 
    # are 100 hyperedges of size 2 and 10 of size 3.
    num_hyperedges::Union{Tuple,Vector{<:Tuple}} = (100, 10)
end

@with_kw mutable struct ModelParams
    # True if the model is solved in discrete time, false if in continuous.
    is_discrete::Union{Bool,Vector{Bool}} = true
    # Which adaptivity rule to use.
    adaptivity_rule::Union{AdaptivityRule,Vector{<:AdaptivityRule}} = RewireToRandom()
    # Which propagation rule to use.
    propagation_rule::Union{PropagationRule,Vector{<:PropagationRule}} = MajorityVoting()
    # Maximum duration of the simulation in continuous time. 
    max_duration::Union{Float64,Vector{Float64}} = 100.0
    # Probability that the adaptivity rule (and not the propagation rule) is executed in a given time step.
    # Only relevant if is_discrete = true
    adaptivity_prob::Union{Float64,Vector{Float64}} = 0.5
    # The rate of the adaptivity events. Only relevant if is_discrete = false.
    adaptivity_rate::Union{Float64,Vector{Float64}} = 1.0
    # The rate of the propagation events. Only relevant if is_discrete = false.
    propagation_rate::Union{Float64,Vector{Float64}} = 1.0
end

@with_kw mutable struct VisualizationParams
    # Take measurements only every n-th time-step. If skip_points = 1, measurements are collected at every time-step. 
    skip_points::Int64 = 1
    # Number of time-steps to simulate before the data is written to file or to the interactive visualization. 
    buffer_size::Int64 = 100
    # The colormap used to plot the states of nodes
    node_colormap::Symbol = :RdYlGn_6
    # The colormap used to plot hyperedges of different sizes
    hyperedge_colormap::Symbol = :thermal
    # The colormap used for any other plots
    misc_colormap::Symbol = :Dark2_3
    # Which panels (or subplots) to display in the interactive visualization. This parameter also determines 
    # which measurements will be collected. Some measurements like the number of motifs are expensive to 
    # collect, so if a measurement is not required by any of the panels, it will not be measured. 
    # Therefore, this parameter is relevant even if the simulation is executed in parallel with MPI 
    # without an interactive visualization. 
    # A map with dependencies of panels on measurements can be found in the file NinjaDashboard.jl.
    # Has to be a symbol with the same name as an AbstractPanel struct. 
    panels::Vector{Symbol} = [:StateDistPanel,
                              :HyperedgeDistPanel,
                              :ActiveHyperedgeDistPanel,
                              :SlowManifoldPanel]
end

@with_kw mutable struct BatchParams
    # If true, a video of the evolution of the dashboard is recorded
    record_video::Bool = false
    # How often to repeat a simulation at the same parameters to collect different stochastic
    # realizations of the trajectory. 
    batch_size::Int64 = 10
    # If true, a prompt will be shown asking if the data should be saved to a file. If false, the data *will not* 
    # be saved if the simulation is execited without MPI and it *will always* be saved if the 
    # simulation is executed in parallel with MPI. 
    prompt_for_save::Bool = false
    # A tag which is appended to the name of the save folder
    save_tag::Union{String,Nothing} = nothing
    # If true, the simulation will be execited non-interactively (without a life visualization) with MPI support. 
    # To actually run the simulation in parallel, you have to jun the main.jl script with mpirun or mpiexec *and* set 
    # with_mpi = true. 
    with_mpi::Bool = false
end

"""
    InputParams

The parameters of the simulation. See the file ./simulation/InputParams.jl for the 
documentation of individual settings. Default values are provided for all settings. 

Every network and model parameter allows to pass a vector of values instead of a single 
value. In this case, simulations are executed for every value in the vector. It is possible 
to vectorize several parameters - in this case, a Cartesian product of all sets is taken. 

# Examples
    InputParams(
        NetworkParams(num_nodes = 42),
        ModelParams(adaptivity_prob = 0.1)
    )

This creates a network with 42 nodes and adaptivity_prob = 0.1. All other parameters 
like the number of hyperedges, initial magnetization, visualization and batch parameters are set to default values. 

    InputParams(
        ModelParams(adaptivity_prob = [0.25, 0.5, 0.75])
    )

This executes three batches of the simulation, one for each value of adaptivity_prob. 
"""
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
Furthermore, only vectors in NetworkParams and ModelParams are supported, since it doesn't make much sense 
to change visualization and batch params between batches.
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

"""
    get_expandable_params(params::InputParams)

Get a list of params which are given in a vector-form and will be expanded in `expand(params)`.
"""
function get_expandable_params(params::InputParams)
    param_names = Dict{Symbol,Vector{Symbol}}()
    param_names[:nparams] = Symbol[]
    param_names[:mparams] = Symbol[]

    nparams, mparams = params.network_params, params.model_params
    for key in fieldnames(NetworkParams)
        if typeof(getfield(nparams, key)) <: Vector
            push!(param_names[:nparams], key)
        end
    end
    for key in fieldnames(ModelParams)
        if typeof(getfield(mparams, key)) <: Vector
            push!(param_names[:mparams], key)
        end
    end
    return param_names
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

function load_params(path::String)
    json_string = read(path, String)
    json = JSON3.read(json_string)

    njson = Dict(json.network_params)
    njson[:num_hyperedges] = Tuple(njson[:num_hyperedges])
    nparams = NetworkParams(; njson...)

    mjson = Dict(json.model_params)
    if :num_time_steps in keys(mjson)
        mjson[:max_duration] = mjson[:num_time_steps] / sum(nparams.num_hyperedges)
        delete!(mjson, :num_time_steps)
    else
        mjson[:max_duration] = Float64(mjson[:max_duration])
    end
    mjson[:adaptivity_rule] = eval(Symbol(mjson[:adaptivity_rule]))()
    mjson[:propagation_rule] = eval(Symbol(mjson[:propagation_rule]))()
    mjson[:adaptivity_rate] = Float64(mjson[:adaptivity_rate])
    mjson[:propagation_rate] = Float64(mjson[:propagation_rate])
    mjson[:adaptivity_prob] = Float64(mjson[:adaptivity_prob])
    mparams = ModelParams(; mjson...)

    vjson = Dict(json.visualization_params)
    for key in [:node_colormap, :hyperedge_colormap, :misc_colormap]
        vjson[key] = Symbol(vjson[key])
    end
    vjson[:panels] = [Symbol(panel) for panel in vjson[:panels]]
    vparams = VisualizationParams(; vjson...)

    bparams = BatchParams(; json.batch_params...)

    params = InputParams(nparams, mparams, vparams, bparams)
    return params
end