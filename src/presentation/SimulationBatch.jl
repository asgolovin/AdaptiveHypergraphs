using JSON3
using StructTypes
using DrWatson
using Dates
using GLMakie
using PCRE

export start_simulation

function start_simulation(params::InputParams)
    nparams, mparams, vparams, bparams = params.network_params,
                                         params.model_params,
                                         params.visualization_params,
                                         params.batch_params
    save_to_file = _prompt_for_save()
    if save_to_file
        output_folder = _create_batch_folder()
        _save_params(params, output_folder)
    end

    n = nparams.num_nodes
    network = HyperNetwork(n, nparams.infected_prob)
    build_RSC_hg!(network, nparams.num_hyperedges)

    model = _create_model(network, mparams)

    dashboard = Dashboard(model; vparams.dashboard_params...)

    for t in 1:(bparams.batch_size)
        reset!(dashboard, model)
        run!(dashboard, mparams.num_time_steps, vparams.steps_per_update)
        sleep(1)
        network = HyperNetwork(n, nparams.infected_prob)
        build_RSC_hg!(network, nparams.num_hyperedges)
        model = _create_model(network, mparams)
    end

    if save_to_file
        rules = "$(typeof(mparams.propagation_rule))_$(typeof(mparams.adaptivity_rule))"
        save(dashboard, output_folder, "$rules.png")
    end
    return nothing
end

function _create_model(network, mparams)
    propagation_rule = mparams.propagation_rule
    adaptivity_rule = mparams.adaptivity_rule
    if mparams.is_discrete
        model = DiscrModel{typeof(propagation_rule),
                           typeof(adaptivity_rule)}(network,
                                                    propagation_rule,
                                                    adaptivity_rule,
                                                    mparams.propagation_prob)
    else
        model = ContinuousModel{typeof(propagation_rule),
                                typeof(adaptivity_rule)}(network,
                                                         propagation_rule,
                                                         adaptivity_rule,
                                                         mparams.propagation_rate,
                                                         mparams.adaptivity_rate)
    end
    return model
end

function _save_params(params::InputParams, folder)
    open(joinpath(folder, "input_params.json"), "w") do io
        return save_json(io, params)
    end
end

function _create_batch_folder()
    println("Enter a tag to name the data folder: ")
    tag = readline()
    filter
    timestamp = Dates.format(now(), "YYYY-mm-dd_HH-MM-SS")
    dir_name = joinpath(projectdir(), "results", "run_$timestamp")
    @assert !ispath(dir_name) "The directory already exists"
    mkpath(dir_name)
    return dir_name
end

function _prompt_for_save()
    println("Save the input parameters and the results of the simulation to a file? [y/n] ")
    s = readline()
    if s == "y"
        dirname = _create_batch_folder()
        println("Ok, saving the input parameters.")
        return true
    else
        println("Ok, the results will not be saved.")
        return false
    end
end