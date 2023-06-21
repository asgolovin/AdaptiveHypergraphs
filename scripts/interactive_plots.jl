include("./interactive_plotting_tools.jl")

# =============================================================
# -------------------------- INPUT ----------------------------

input_folder = joinpath(projectdir(),
                        "data/run_2023-04-24_14-46-47_maj_voting_rts")

panel_symbols = [:StateDistPanel, :ActiveHyperedgeDistPanel, :SlowManifoldPanel]

skip_points = 10

# -------------------------------------------------------------
# =============================================================

fig, grid_pos = create_fig(panel_symbols)
panels, measurements = create_panels(joinpath(input_folder, "batch_001", "run_001"),
                                     panel_symbols,
                                     grid_pos)

for batchdir in readdir(input_folder; join=true)
    if !isdir(batchdir)
        continue
    end
    batch_num = parse(Int64, match(r"batch_([0-9]+)", batchdir)[1])
    if batch_num != 17
        continue
    end

    input_params = load_params(joinpath(batchdir, "input_params.json"))
    @show input_params.model_params.adaptivity_prob

    for rundir in readdir(batchdir; join=true)
        if !isdir(rundir)
            continue
        end

        run_num = parse(Int64, match(r"run_([0-9]+)", splitdir(rundir)[end])[1])

        for panel in panels
            deactivate_lines!(panel)
        end

        for meas_sym in keys(measurements)
            load_meas(rundir, measurements[meas_sym])
        end

        for meas_sym in keys(measurements)
            for meas in measurements[meas_sym]
                notify(meas.log.observable_values)
            end
        end
    end
end

display(fig)