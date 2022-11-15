using JSON3
using StructTypes
using DrWatson
using Dates
using REPL.TerminalMenus

export start_simulation

function start_simulation(params::InputParams)
    bparams = params.batch_params
    vparams = params.visualization_params

    # turns on a prompt if the data should be saved. 
    if bparams.prompt_for_save
        save_to_file = _prompt_for_save()
        if save_to_file
            output_folder = _create_batch_folder()
            _save_params(params, output_folder)
        end
    else
        save_to_file = false
    end

    # a vector with all possible combinations of all sweeped params
    param_vector = expand(params)

    nparams = param_vector[1].network_params
    mparams = param_vector[1].model_params

    n = nparams.num_nodes
    network = HyperNetwork(n, nparams.infected_prob)
    build_RSC_hg!(network, nparams.num_hyperedges)

    model = _create_model(network, mparams)

    dashboard = Dashboard(model; vparams)

    for param in param_vector
        nparams = param.network_params
        mparams = param.model_params

        for t in 1:(bparams.batch_size)
            if t != 1
                network = HyperNetwork(n, nparams.infected_prob)
                build_RSC_hg!(network, nparams.num_hyperedges)
                model = _create_model(network, mparams)
                reset!(dashboard, model)
            end
            run!(dashboard, mparams.num_time_steps)
            sleep(0.1)
        end
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
                                                    mparams.adaptivity_prob)
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
    # replace all spaces by underscores and throw out all non-alphanumeric characters
    tag = replace(tag, " " => "_", r"[^\p{L}\p{N},_]" => "")
    timestamp = Dates.format(now(), "YYYY-mm-dd_HH-MM-SS")
    dir_name = joinpath(projectdir(), "results", "run_$(timestamp)_$(tag)")
    @assert !ispath(dir_name) "The directory already exists"
    mkpath(dir_name)
    local_dir = joinpath("..", "results", "run_$(timestamp)_$(tag)")
    println("""The results will be saved to the folder "$local_dir".""")
    return dir_name
end

function _prompt_for_save()
    ans = Base.prompt("Press space and then twice Enter")
    options = ["no", "yes"]
    menu = RadioMenu(options)
    choice = request("Save the input parameters and the results of the simulation to a file?",
                     menu)
    if choice == 2
        return true
    else
        println("Ok, the results will not be saved.")
        return false
    end
end