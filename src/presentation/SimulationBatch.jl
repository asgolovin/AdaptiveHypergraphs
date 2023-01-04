using JSON3
using StructTypes
using DrWatson
using Dates
using REPL.TerminalMenus
using MPI

export start_simulation

function start_simulation(params::InputParams)
    bparams = params.batch_params
    vparams = params.visualization_params

    if bparams.with_mpi
        MPI.Init()
        comm = MPI.COMM_WORLD
        rank = MPI.Comm_rank(comm)
        size = MPI.Comm_size(comm)
    else
        rank = 0
        size = 1
    end

    # turns on a prompt if the data should be saved. 
    save_to_file = false
    if bparams.prompt_for_save
        if rank == 0
            save_to_file = _prompt_for_save()
            if save_to_file
                output_folder = _create_batch_folder()
                _save_params(params, output_folder)
            end
        end

        # [MPI] broadcast the answer to all ranks
        if bparams.with_mpi
            save_to_file = MPI.Bcast!(save_to_file, 0, comm)
        end
    else
        save_to_file = false
    end

    # a vector of params which can be expanded
    expandable_params = get_expandable_params(params)
    # a vector with all possible combinations of all sweeped params
    param_vector = expand(params)

    nparams = param_vector[1].network_params
    mparams = param_vector[1].model_params

    n = nparams.num_nodes
    max_size = length(nparams.num_hyperedges) + 1
    network = HyperNetwork(n, nparams.infected_prob, max_size)
    build_RSC_hg!(network, nparams.num_hyperedges)

    model = _create_model(network, mparams)

    # create an invisible dashboard on all ranks with MPI and a normal one without
    if bparams.with_mpi
        dashboard = NinjaDashboard(model, vparams)
    else
        dashboard = Dashboard(model; vparams)
    end

    for (i, param) in enumerate(param_vector)
        nparams = param.network_params
        mparams = param.model_params

        if rank == 0
            println("\nExecuting batch $i/$(length(param_vector))")
            for param in expandable_params[:nparams]
                println("$param = $(getfield(nparams, param))")
            end
            for param in expandable_params[:mparams]
                println("$param = $(getfield(mparams, param))")
            end
        end

        for t in 1:(bparams.batch_size)
            # Split the work among MPI ranks
            if mod((i - 1) * bparams.batch_size + t, size) != rank
                continue
            end

            num_batches = length(param_vector)
            batch_size = bparams.batch_size
            println("[rank $rank] executing batch $i/$num_batches, iteration $t/$batch_size")

            max_size = length(nparams.num_hyperedges) + 1
            network = HyperNetwork(n, nparams.infected_prob, max_size)
            build_RSC_hg!(network, nparams.num_hyperedges)
            model = _create_model(network, mparams)
            reset!(dashboard, model)

            run!(dashboard, mparams.num_time_steps)
            sleep(0.01)
        end
    end

    if save_to_file
        save(dashboard, output_folder, "dashboard.png")
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