include("figure_plotting_tools.jl")

# ======================================================
# -------------------- INPUT ---------------------------

folder_maj_rts = joinpath(projectdir(),
                          "final_data/run_2023-02-06_14-54-33_maj_voting_rts")
folder_maj_rtr = joinpath(projectdir(),
                          "final_data/run_2023-02-06_14-54-34_maj_voting_rtr")
folder_prop_rts = joinpath(projectdir(),
                           "final_data/run_2023-02-07_00-09-40_prop_voting_rts_merged")
folder_prop_rtr = joinpath(projectdir(),
                           "final_data/run_2023-02-06_21-50-23_prop_voting_rtr")

panels = [:maj_rts, :maj_rtr, :prop_rts, :prop_rtr]
titles = Dict(:maj_rts => "Majority voting, rewire-to-same",
              :maj_rtr => "Majority voting, rewire-to-random",
              :prop_rts => "Prop. voting, rewire-to-same",
              :prop_rtr => "Prop. voting, rewire-to-random")
input_folders = Dict(:maj_rts => folder_maj_rts,
                     :maj_rtr => folder_maj_rtr,
                     :prop_rts => folder_prop_rts,
                     :prop_rtr => folder_prop_rtr)

num_nodes = 10000
max_size = 4
hyperedge_colormap = :thermal
linecolors = get(colorschemes[hyperedge_colormap], [1, 2.5, 3], (1, max_size))
# number of datapoints that are skipped for the fits
skip_initial = Dict(:maj_rts => 1,
                    :maj_rtr => 1,
                    :prop_rts => 300,
                    :prop_rtr => 100)
# plot every nth point 
skip_middle = Dict(:maj_rts => 1,
                   :maj_rtr => 1,
                   :prop_rts => 5,
                   :prop_rtr => 10)
plot_fits = false
prompt = true

# which run to plot in bold
highlight_run = Dict(:maj_rts => [5, 5, 5, 5, 5],
                     :maj_rtr => [5, 5, 5, 5, 5],
                     :prop_rts => [9, 5, 3, 3, 7],
                     :prop_rtr => [5, 5, 5, 5, 5])

filename = "./figures/thesis/slow_manifold_hypergraph_all_four.pdf"

# ------------------------------------------------------
# ======================================================

# create figure
size_inches = (6, 4)
size_pt = 1.5 * 72 .* size_inches
fig = Figure(; resolution=size_pt)

axis = Dict()
for (i, panel) in enumerate(panels)
    row = (i + 1) ÷ 2
    col = mod1(i, 2)
    ax_gridpos = fig[1, 1]
    axis[panel] = Axis(ax_gridpos[row, col];
                       titlesize=15,
                       title=titles[panel],
                       xticklabelsize=12,
                       yticklabelsize=12,
                       xlabelsize=13,
                       ylabelsize=13,
                       titlefont="Latin Modern Roman")
    if row == 1
        hidexdecorations!(axis[panel]; grid=false)
    else
        axis[panel].xlabel = L"magnetization $m$"
    end
    if col == 2
        hideydecorations!(axis[panel]; grid=false)
    else
        axis[panel].ylabel = L"density of active links $\rho_{d_i}$"
    end
    axis[panel].yticks = 0.0:0.2:0.6
    #axis[panel].xminorticks = -1:0.1:1
    #axis[panel].xminorticksvisible = true
    #axis[panel].xminorgridvisible = true
end

for panel in panels
    magnetization_data = []
    density_data = Dict(s => [] for s in 2:max_size)

    data_folder = DataFolder(input_folders[panel], :simulation)
    ax = axis[panel]
    skip_init = skip_initial[panel]
    skip_mid = skip_middle[panel]

    for (batch_folder, batch_num) in data_folder
        for (run_folder, run_num) in batch_folder
            if panel == :prop_rts && run_num > 5
                skip_mid = 100
            end

            filename = "state_count_A.csv"
            path = joinpath(run_folder, filename)
            df = CSV.read(path, DataFrame; delim=", ")

            state_count_A = df.value
            magnetization = 2 .* (state_count_A ./ num_nodes) .- 1
            magnetization = vcat(magnetization[1:skip_init],
                                 magnetization[(skip_init + 1):skip_mid:end])

            if batch_num == 3
                append!(magnetization_data, magnetization[skip_init:end])
            end

            for size in 2:max_size
                filename = "active_hyperedge_count_$size.csv"
                path = joinpath(run_folder, filename)
                df = CSV.read(path, DataFrame; delim=", ")

                active_hyperedeges = df.value
                active_density = active_hyperedeges ./ num_nodes
                active_density = vcat(active_density[1:skip_init],
                                      active_density[(skip_init + 1):skip_mid:end])

                if batch_num == 3
                    append!(density_data[size], active_density[skip_init:end])
                end

                linecolor = linecolors[size - 1]
                label = "size $size"

                # plot one trajectory in a darker color and with a label
                if run_num == highlight_run[panel][batch_num]
                    lines!(ax, magnetization, active_density; color=linecolor,
                           label=label, linewidth=0.7)
                else
                    lines!(ax, magnetization, active_density; color=(linecolor, 0.3),
                           linewidth=0.7)
                end
            end
        end
    end

    xlims!(ax, (-1, 1))
    ylims!(ax, (0, 0.65))

    if plot_fits
        if panel != :maj_rts && panel != :maj_rtr
            for size in 2:max_size
                p, cov = fit_parabola(magnetization_data, density_data[size], false)
                x = collect(-1:0.01:1)
                lines!(ax, x, p.(x); color=:red, label="polynomial fit")
            end
        end
    end
end

leg = Legend(fig[1, 2], axis[:prop_rtr]; merge=true, orientation=:vertical, labelsize=12)

display(fig)

if prompt
    prompt_for_save(filename, fig; pt_per_unit=0.666)
end