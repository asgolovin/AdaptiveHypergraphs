include("figure_plotting_tools.jl")

# ======================================================
# -------------------- INPUT ---------------------------

input_folder = joinpath(projectdir(),
                        "final_data/run_2023-02-25_16-12-33_maj_voting_rts_p01")

panels = [:maj_rts, :maj_rtr, :prop_rts, :prop_rtr]
title = "Majority voting, rewire-to-same"

num_nodes = 10000
max_size = 4
hyperedge_colormap = :thermal
linecolors = get(colorschemes[hyperedge_colormap], [1, 2.5, 3], (1, max_size))

plot_fits = false
prompt = true

# which run to plot in bold
highlight_run = [5, 5, 1, 5, 5]

filename = "./figures/slow_manifold_hypergraph_special_p01.pdf"

# ------------------------------------------------------
# ======================================================

# create figure
size_inches = (3.7, 2.3)
size_pt = 1.5 * 72 .* size_inches
fig = Figure(; resolution=size_pt)

ax = Axis(fig[1, 1];
          titlesize=15,
          title=title,
          xticklabelsize=12,
          yticklabelsize=12,
          xlabelsize=13,
          ylabelsize=13,
          titlefont="Latin Modern Roman")

ax.xlabel = L"magnetization $m$"

ax.ylabel = L"density of active links $\rho_{d_i}$"
ax.yticks = 0.0:0.2:0.6
#ax[panel].xminorticks = -1:0.1:1
#ax[panel].xminorticksvisible = true
#ax[panel].xminorgridvisible = true

magnetization_data = []
density_data = Dict(s => [] for s in 2:max_size)

data_folder = DataFolder(input_folder, :simulation)

for (batch_folder, batch_num) in data_folder
    for (run_folder, run_num) in batch_folder
        filename = "state_count_A.csv"
        path = joinpath(run_folder, filename)
        df = CSV.read(path, DataFrame; delim=", ")

        state_count_A = df.value
        magnetization = 2 .* (state_count_A ./ num_nodes) .- 1

        for size in 2:max_size
            filename = "active_hyperedge_count_$size.csv"
            path = joinpath(run_folder, filename)
            df = CSV.read(path, DataFrame; delim=", ")

            active_hyperedeges = df.value
            active_density = active_hyperedeges ./ num_nodes

            linecolor = linecolors[size - 1]
            label = "size $size"

            # plot one trajectory in a darker color and with a label
            if run_num == highlight_run[batch_num]
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

leg = Legend(fig[1, 2], ax; merge=true, orientation=:vertical, labelsize=12)

display(fig)

if prompt
    prompt_for_save(filename, fig; pt_per_unit=0.666)
end