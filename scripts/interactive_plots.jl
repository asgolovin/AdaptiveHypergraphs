include("./interactive_plotting_tools.jl")

# =============================================================
# -------------------------- INPUT ----------------------------

input_folder = joinpath(projectdir(),
                        "final_data/run_2023-02-23_13-32-26_motifs_D4")

panel_symbols = [:StateDistPanel, :HyperedgeDistPanel, :MomentClosurePanel]

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

    for rundir in readdir(batchdir; join=true)
        if !isdir(rundir)
            continue
        end

        run_num = parse(Int64, match(r"run_([0-9]+)", splitdir(rundir)[end])[1])

        if run_num != 1
            continue
        end

        println("$run_num")

        for panel in panels
            deactivate_lines!(panel)
        end

        for meas_sym in keys(measurements)
            println("$meas_sym")
            load_meas(rundir, measurements[meas_sym])
        end

        for meas_sym in keys(measurements)
            for meas in measurements[meas_sym]
                println("$meas_sym")
                notify(meas.log.observable_values)
            end
        end
    end
end

display(fig)