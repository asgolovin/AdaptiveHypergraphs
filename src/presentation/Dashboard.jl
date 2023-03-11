export Dashboard, run!, record!, reset!

"""
    Dashboard <: AbstractDashboard

A dashboard that actually visualizes stuff.
"""
struct Dashboard <: AbstractDashboard
    fig::Figure
    panels::Vector{AbstractPanel}
    mo::ModelObservable
    measurement_types::Vector{DataType}
end

#! format: off
"""
    Dashboard(model::AbstractModel;
              vparams::VisualizationParams)

A visualization of the evolution of the hypergraph during the simulation.

# Arguments
- `model::AbstractModel` - the model on which the dashboard is based on. 
- `vparams::VisualizationParams` - parameters used for visualization
"""
function Dashboard(model::AbstractModel,
                   vparams::VisualizationParams; save_folder::Union{Nothing,String})
    #! format: on
    fig = Figure(; resolution=(1200, 800))
    display(fig)

    # collect a list of the measurements on which the panels depend on
    panel_symbols = vparams.panels
    measurements = Vector{DataType}()
    for panel in panel_symbols
        append!(measurements, PANEL_DEPENDENCIES[panel])
    end
    measurement_types = unique(measurements)

    mo = ModelObservable(model, measurement_types;
                         skip_points=vparams.skip_points,
                         buffer_size=vparams.buffer_size,
                         write_to_observables=true,
                         save_folder=save_folder)

    # Display plots in the left part of the figure and the info box on the right. 
    plot_box = fig[1, 1] = GridLayout()
    # info_box = fig[1, 2] = GridLayout()

    # ==========================================================================================
    # -------------------------------------- PLOTS ---------------------------------------------

    # determine the number of rows and columns 
    num_panels = length(panel_symbols)
    nrows = Int64(floor(sqrt(num_panels)))
    ncols = Int64(ceil(sqrt(num_panels)))

    panels = []
    graph_properties = Dict(:num_nodes => get_num_nodes(mo.network[]),
                            :max_hyperedge_size => get_max_size(mo.network[]),
                            :num_hyperedges => get_num_hyperedges(mo.network[]))

    for (i, panel) in enumerate(panel_symbols)
        panel_type = eval(panel)
        col = mod1(i, ncols)
        row = (i - 1) ÷ ncols + 1

        # HypergraphPanel doesn't take a Measurement and instead needs the network, so 
        # we treat it separately
        if panel_type <: HypergraphPanel
            panel = HypergraphPanel(plot_box[col, row], mo.network,
                                    graph_properties,
                                    vparams)
        else
            # add the dependent measurements to the arg list
            measurement_types = PANEL_DEPENDENCIES[panel]
            measurements = Dict()
            for type in measurement_types
                sym = Symbol(_snake_case("$type"))
                measurements[sym] = getfield(mo, sym)
            end
            panel = panel_type(plot_box[col, row], measurements, graph_properties, vparams)
        end
        push!(panels, panel)
    end

    # ==========================================================================================
    # ------------------------------------ INFO BOX---------------------------------------------

    # TODO

    return Dashboard(fig, panels, mo, measurement_types)
end

"""
    record!(dashboard::Dashboard, filename::String, num_steps::Int64, framerate::Int64)

Run the simulation for `num_steps` time steps and record a video of the dashboard. 

The video is saved to ./videos in a .mp4 format. `filename` should only contain the name 
of the file without the extension.
"""
function record!(dashboard::Dashboard, filename::String, num_steps::Int64,
                 framerate::Int64)
    savepath = joinpath("videos", filename * ".mp4")

    num_updates = num_steps ÷ mo.buffer_size
    record(dashboard.fig, savepath, 1:num_updates; framerate=framerate, compression=1) do i
        return run!(dashboard, mo.buffer_size)
    end
    return dashboard
end

"""
    reset!(dashboard::Dashboard, model::AbstractModel)

Reset the dashboard to run the next simulation from the batch. 

The old history plot lines are made inactive and are grayed out. 
The observables in the ModelObservable are reset to track the new data from `model`.
"""
function reset!(dashboard::AbstractDashboard, model::AbstractModel,
                save_folder::Union{Nothing,String})
    # gray out the history plot lines
    for panel in dashboard.panels
        if typeof(panel) <: AbstractTimeSeriesPanel ||
           typeof(panel) <: SlowManifoldPanel
            deactivate_lines!(panel)
        end
    end

    # reset observables
    rebind_model!(dashboard.mo, model, save_folder)
    return dashboard
end

function set_solution(dashboard::Dashboard, t::Vector{Float64},
                      sol::Dict{AbstractMotif,Vector{Float64}})
    for panel in dashboard.panels
        set_solution(panel, t, sol)
    end
end

"""
    save(dashboard::Dashboard, folder::String, filename::String)

Save the state of the dashboard as a figure.
The filename should be given with an extension, for example, dash.png. 
"""
function save(dashboard::Dashboard, folder::String, filename::String)
    return GLMakie.save(joinpath(folder, filename), dashboard.fig)
end