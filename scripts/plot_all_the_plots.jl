include("./plotting_tools.jl")

# =============================================================
# -------------------------- INPUT ----------------------------

input_folder = joinpath(projectdir(),
                        "results/run_2023-02-06_14-54-33_maj_voting_rts")

panel_symbols = [:StateDistPanel, :HyperedgeDistPanel, :ActiveHyperedgeDistPanel,
                 :SlowManifoldPanel]

skip_points = 100

# -------------------------------------------------------------
# =============================================================

grid_pos = create_fig(panel_symbols)
panels, measurements = create_panels(joinpath(input_folder, "batch_001", "run_001"),
                                     panel_symbols,
                                     grid_pos)

for batchdir in readdir(input_folder; join=true)
    if !isdir(batchdir)
        continue
    end
    batch_num = parse(Int64, match(r"batch_([0-9]+)", batchdir)[1])

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
    end
end