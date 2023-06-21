using CairoMakie
using DrWatson
using CSV
using Colors
using ColorSchemes
using DataFrames

using AdaptiveHypergraphs

GLMakie.activate!()

Label = AdaptiveHypergraphs.Label

# structure of the folder: 
# input_folder
# | batch_001
# | | input_params.json
# | | run_001
# | | | <measurement_name>_<label>.csv
# | | | ...
# | | run_002
# | | <run_measurement_name>_<label>.csv

hyperedge_colormap = :thermal
function hyperedge_linecolors(max_size)
    return get(colorschemes[hyperedge_colormap], [1, 2.5, 3],
               (1, max_size))
end

"""
Create a vector of empty measurements of a specific type
"""
function create_measurements(folder, meas_type::DataType)
    meas_str = AdaptiveHypergraphs._snake_case("$meas_type")

    # get a list of labels from file names
    filenames = readdir(folder)
    filter!(filename -> startswith(filename, meas_str), filenames)
    files = [joinpath(folder, name) for name in filenames]

    measurements = meas_type[]

    for file in files
        m = match(r"([a-z_]+)_([a-zA-Z0-9\[\]\|\h]+).csv", file)
        if isnothing(m)
            throw(Exception("The name of the file does not match the required pattern."))
        end

        label_str = String(m[2])
        labeltype = fieldtype(meas_type, :label)
        if labeltype <: State
            label = label_str == "A" ? AdaptiveHypergraphs.A : AdaptiveHypergraphs.B
        elseif labeltype <: AbstractMotif
            label = motif_from_string(label_str)
        elseif labeltype <: Int64
            label = parse(Int64, label_str)
        else
            println("Encountered an unknown label type $labeltype. Treating it as a string.")
            label = label_str
        end

        meas = meas_type(label; save_folder=nothing)

        push!(measurements, meas)
    end
    return measurements
end

"""
Load data into a vector of measurements
"""
function load_meas(folder, measurements::Vector)
    meas_type = typeof(measurements[1])
    meas_str = AdaptiveHypergraphs._snake_case("$meas_type")

    for meas in measurements
        label = meas.label
        filename = "$(meas_str)_$label.csv"
        path = joinpath(folder, filename)
        df = CSV.read(path, DataFrame; delim=", ")
        empty!(meas.log.observable_indices[])
        empty!(meas.log.observable_values[])
        append!(meas.log.observable_indices[], df.index[1:skip_points:end])
        append!(meas.log.observable_values[], df.value[1:skip_points:end])
    end

    return measurements
end

function create_fig(panel_symbols)
    size_inches = (10, 6)
    size_pt = 72 .* size_inches
    fig = Figure(; resolution=size_pt, fontsize=16)

    num_axes = length(panel_symbols)
    nrows = Int64(floor(sqrt(num_axes)))
    ncols = Int64(ceil(sqrt(num_axes)))
    grid = Dict()

    for (i, panel_sym) in enumerate(panel_symbols)
        col = mod1(i, ncols)
        row = (i - 1) ÷ ncols + 1
        grid[panel_sym] = fig[row, col]
    end
    return fig, grid
end

"""
Create panels in the figure
"""
function create_panels(folder, panel_symbols, grid)
    panels = []
    batchdir = joinpath(splitpath(folder)[1:(end - 1)])
    params = load_params(joinpath(batchdir, "input_params.json"))

    measurements = Dict()
    for panel_sym in panel_symbols
        # filter out the measurements needed for a panel
        meas_types = PANEL_DEPENDENCIES[panel_sym]
        for meas_type in meas_types
            if !(meas_type <: AbstractStepMeasurement)
                continue
            end
            meas_sym = Symbol(AdaptiveHypergraphs._snake_case("$meas_type"))
            if meas_sym in keys(measurements)
                continue
            end
            # create a vector of Measurements from the data
            meas_vector = create_measurements(folder, meas_type)
            measurements[meas_sym] = meas_vector
        end

        # give the Measurements to the functions from Panels.jl
        panel_type = eval(panel_sym)

        nparams = params.network_params
        max_size = length(nparams.num_hyperedges) + 1
        num_hyperedges = sum(nparams.num_hyperedges)
        graph_properties = Dict(:num_nodes => nparams.num_nodes,
                                :max_hyperedge_size => max_size,
                                :num_hyperedges => num_hyperedges)
        panel = panel_type(grid[panel_sym],
                           measurements,
                           graph_properties,
                           params.visualization_params)
        push!(panels, panel)
    end
    return panels, measurements
end