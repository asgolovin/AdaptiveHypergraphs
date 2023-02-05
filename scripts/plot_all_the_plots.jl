using GLMakie
using DrWatson
using CSV
using DataFrames

using AdaptiveHypergraphs

Label = AdaptiveHypergraphs.Label

# =============================================================
# -------------------------- INPUT ----------------------------

input_folder = joinpath(projectdir(),
                        "results/run_2023-02-04_13-02-16_maj_voting_rtr")

panels = [:StateDistPanel, :HyperedgeDistPanel, :ActiveHyperedgeDistPanel,
          :SlowManifoldPanel, :MomentClosurePanel]

skip_points = 100

# -------------------------------------------------------------
# =============================================================

# structure of the folder: 
# input_folder
# | batch_001
# | | input_params.json
# | | run_001
# | | | <measurement_name>_<label>.csv
# | | | ...
# | | run_002
# | | <run_measurement_name>_<label>.csv

"""
Create a vector of Measurements from a given folder with data
"""
function load_meas(folder, meas_type::DataType)
    filenames = readdir(folder)
    meas_str = AdaptiveHypergraphs._snake_case("$meas_type")
    filter!(filename -> startswith(filename, meas_str), filenames)
    files = [joinpath(folder, name) for name in filenames]

    measurements = meas_type[]

    for file in files
        m = match(r"([a-z_]+)_([a-zA-Z0-9\[\]\|\h]+).csv", file)
        if isnothing(m)
            throw(Exception("The name of the file does not match the required pattern."))
        end
        labeltype = meas_type
        label_str = String(m[2])
        labeltype = fieldtype(meas_type, :label)
        if labeltype <: State
            label = label_str == "A" ? AdaptiveHypergraphs.A : AdaptiveHypergraphs.B
        elseif labeltype <: Label
            label = Label(label_str)
        elseif labeltype <: Int64
            label = parse(Int64, label_str)
        else
            println("Encountered an unknown label type $labeltype. Treaning it as a string.")
            label = label_str
        end

        df = CSV.read(file, DataFrame; delim=", ")
        IndexType = typeof(df.index[1])
        ValueType = typeof(df.value[1])

        indices = df.index[1:100:end]
        values = df.value[1:100:end]
        log = MeasurementLog{IndexType,ValueType}(Observable(indices),
                                                  Observable(values))
        push!(measurements, meas_type(log, label))
    end

    return measurements
end

for batchdir in readdir(input_folder; join=true)
    if !isdir(batchdir)
        continue
    end
    batch_num = parse(Int64, match(r"batch_([0-9]+)", batchdir)[1])
    fig = Figure()
    display(fig)

    num_axes = length(panels)
    nrows = Int64(floor(sqrt(num_axes)))
    ncols = Int64(ceil(sqrt(num_axes)))
    grid_pos = Dict()

    for (i, panel) in enumerate(panels)
        col = mod1(i, ncols)
        row = (i - 1) ÷ ncols + 1
        grid_pos[panel] = fig[row, col]
    end

    # get input params
    params = load_params(joinpath(batchdir, "input_params.json"))

    for rundir in readdir(batchdir; join=true)
        if !isdir(rundir)
            continue
        end

        run_num = parse(Int64, match(r"run_([0-9]+)", splitdir(rundir)[end])[1])

        println(run_num)

        #if run_num != 1
        #    continue
        #end

        files = readdir(rundir)

        for panel_sym in panels
            # filter out the measurements needed for a panel
            meas_types = PANEL_DEPENDENCIES[panel_sym]
            measurements = Dict()
            for meas_type in meas_types
                # create a vector of Measurements from the data
                meas_vector = load_meas(rundir, meas_type)
                measurements[Symbol(AdaptiveHypergraphs._snake_case("$meas_type"))] = meas_vector
            end

            # give the Measurements to the functions from Panels.jl
            panel_type = eval(panel_sym)

            nparams = params.network_params
            max_size = length(nparams.num_hyperedges) + 1
            num_hyperedges = sum(nparams.num_hyperedges)
            graph_properties = Dict(:num_nodes => nparams.num_nodes,
                                    :max_hyperedge_size => max_size,
                                    :num_hyperedges => num_hyperedges)
            panel = panel_type(grid_pos[panel_sym],
                               measurements,
                               graph_properties,
                               params.visualization_params)
            deactivate_lines!(panel)
        end
    end
end