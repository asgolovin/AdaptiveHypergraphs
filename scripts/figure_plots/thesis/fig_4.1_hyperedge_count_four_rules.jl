include("figure_plotting_tools.jl")

# ======================================================
# -------------------- INPUT ---------------------------

folder_maj_rts = joinpath(projectdir(),
                          "data/thesis/run_2023-02-06_14-54-33_maj_voting_rts")
folder_maj_rtr = joinpath(projectdir(),
                          "data/thesis/run_2023-02-06_14-54-34_maj_voting_rtr")
folder_prop_rts = joinpath(projectdir(),
                           "data/thesis/run_2023-02-07_00-09-40_prop_voting_rts_merged")
folder_prop_rtr = joinpath(projectdir(),
                           "data/thesis/run_2023-02-06_21-50-23_prop_voting_rtr")

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

# plot every nth point 
skip_middle = Dict(:maj_rts => 1,
                   :maj_rtr => 1,
                   :prop_rts => 5,
                   :prop_rtr => 10)

prompt = true

# which run to plot in bold
highlight_run = Dict(:maj_rts => 5,
                     :maj_rtr => 5,
                     :prop_rts => 5,
                     :prop_rtr => 5)

filename = "./figures/thesis/hyperedge_dist_all_four.pdf"

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
    axis[panel].xticks = 0:20:100
    axis[panel].yticks = 0:2500:12000

    if row == 1
        hidexdecorations!(axis[panel]; grid=false)
    else
        axis[panel].xlabel = L"time $t$"
    end

    if col == 2
        hideydecorations!(axis[panel]; grid=false)
    else
        axis[panel].ylabel = "# of hyperedges"
    end
end

for panel in panels
    data_folder = DataFolder(input_folders[panel], :simulation)
    ax = axis[panel]
    skip_mid = skip_middle[panel]

    for (batch_folder, batch_num) in data_folder
        if batch_num != 3
            continue
        end
        for (run_folder, run_num) in batch_folder
            for size in 2:max_size
                filename = "hyperedge_count_$size.csv"
                path = joinpath(run_folder, filename)
                df = CSV.read(path, DataFrame; delim=", ")

                time = df.index
                hyperedge_count = df.value

                linecolor = linecolors[size - 1]
                label = "size $size"

                # plot one trajectory in a darker color and with a label
                if run_num == highlight_run[panel]
                    lines!(ax, time, hyperedge_count; color=linecolor,
                           label=label, linewidth=1.5)
                else
                    lines!(ax, time, hyperedge_count; color=(linecolor, 0.3),
                           linewidth=1.5)
                end
            end
        end
    end

    xlims!(ax, (0, 100))
    ylims!(ax, (0, 12400))
end

leg = Legend(fig[1, 2], axis[:prop_rtr]; merge=true, orientation=:vertical, labelsize=12)

display(fig)

if prompt
    prompt_for_save(filename, fig; pt_per_unit=0.666)
end