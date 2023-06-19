include("figure_plotting_tools.jl")

# ======================================================
# -------------------- INPUT ---------------------------

folder_maj_rts = joinpath(projectdir(),
                          "data/run_2023-05-03_16-54-25_slow_manifold_maj_voting_rts")
folder_maj_rtr = joinpath(projectdir(),
                          "data/run_2023-05-03_16-54-25_slow_manifold_maj_voting_rtr")
folder_prop_rts = joinpath(projectdir(),
                           "data/run_2023-05-04_08-11-45_slow_manifold_prop_voting_rts")
folder_prop_rtr = joinpath(projectdir(),
                           "data/run_2023-05-03_21-37-37_slow_manifold_prop_voting_rtr")

panels = [:maj_rts, :maj_rtr, :prop_rts, :prop_rtr]
titles = Dict(:maj_rts => "Majority voting,\nrewire-to-same",
              :maj_rtr => "Majority voting,\nrewire-to-random",
              :prop_rts => "Proportional voting,\nrewire-to-same",
              :prop_rtr => "Proportional voting,\nrewire-to-random")
input_folders = Dict(:maj_rts => folder_maj_rts,
                     :maj_rtr => folder_maj_rtr,
                     :prop_rts => folder_prop_rts,
                     :prop_rtr => folder_prop_rtr)

num_nodes = 10000
max_size = 4
linecolors = hyperedge_linecolors(max_size)

# number of datapoints after which the simulation has converged to a parabola
skip_initial = Dict(:maj_rts => 50,
                    :maj_rtr => 50,
                    :prop_rts => 300,
                    :prop_rtr => 200)
# plot every n-th point to reduce the size of the figures
skip_middle = Dict(:maj_rts => 10,
                   :maj_rtr => 10,
                   :prop_rts => 300,
                   :prop_rtr => 300)
plot_fits = false
prompt = true

# which run to plot in bold
highlight_run = Dict(:maj_rts => [5, 5, 5, 5, 5],
                     :maj_rtr => [5, 5, 5, 5, 5],
                     :prop_rts => [9, 5, 3, 3, 7],
                     :prop_rtr => [5, 5, 5, 5, 5])

xlabel = L"m"
ylabel = L"\rho_i"
xticks = -1:0.5:1
yticks = 0:0.2:0.6

filename = "./figures/fig_3_slow_manifold.pdf"

# ------------------------------------------------------
# ======================================================

# create figure
fig = create_figure(:large, 3 / 2)
ax_gridpos = fig[1, 1] = GridLayout()

axis = add_four_rules_axes!(fig, titles, xlabel, ylabel, xticks, yticks)

for panel in panels
    magnetization_data = []
    density_data = Dict(s => [] for s in 2:max_size)

    data_folder = DataFolder(input_folders[panel], :simulation)
    local ax = axis[panel]
    skip_init = skip_initial[panel]
    skip_mid = skip_middle[panel]

    for (batch_folder, batch_num) in data_folder
        for (run_folder, run_num) in batch_folder
            local filename = "state_count_A.csv"
            path = joinpath(run_folder, filename)
            _, state_count_A = load_data(path, run_folder)

            local magnetization = 2 .* (state_count_A ./ num_nodes) .- 1
            magnetization = vcat(magnetization[1:skip_init],
                                 magnetization[(skip_init + 1):skip_mid:end])

            if batch_num == 3
                append!(magnetization_data, magnetization[skip_init:end])
            end

            for size in 2:max_size
                local filename = "active_hyperedge_count_$size.csv"
                path = joinpath(run_folder, filename)
                _, active_hyperedeges = load_data(path, run_folder)

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

# Draw the scetch in the last column

ax_graphical = Dict(:maj => Axis(ax_gridpos[1, 3]),
                    :prop => Axis(ax_gridpos[2, 3]))

ax_graphical[:maj].title = "Majority voting\nphase space sketch"
ax_graphical[:prop].title = "Proportional voting\nphase space sketch"

for ax in values(ax_graphical)
    hidexdecorations!(ax)
    hideydecorations!(ax)
    #xlims!(ax, (-1, 1))
    #ylims!(ax, (0, 2.1))
end

x = collect(-1:0.1:1)
parabola_sketch = -(x .- 1) .* (x .+ 1)
line_kwargs = Dict(:color => :black, :linewidth => 1)
scatter_kwargs = Dict(:color => :black, :marker => :rtriangle)

# draw the parabola and the vertical lines
for ax in values(ax_graphical)
    lines!(ax, x, parabola_sketch; line_kwargs...)
    lines!(ax, x, parabola_sketch; line_kwargs...)
    lines!(ax, [0, 0], [2, 1]; line_kwargs...)
    lines!(ax, [0, 0], [0.1, 1]; line_kwargs...)
    lines!(ax, [-1, 1], [0, 0]; line_kwargs...)
    scatter!(ax, [0], [1.5]; rotations=-π / 2, scatter_kwargs...)
    scatter!(ax, [0], [0.5]; rotations=π / 2, scatter_kwargs...)
    scatter!(ax, [-1, 1], [0, 0]; marker=:circle,
             color=:black, markersize=10)
end

# MAJORITY VOTING
# arrows on the parabola
ind = length(x) - 6
angle = atan(parabola_sketch[ind] - parabola_sketch[ind - 1], x[ind] - x[ind - 1])
scatter!(ax_graphical[:maj], [x[ind], -x[ind]],
         [parabola_sketch[ind], parabola_sketch[ind]];
         rotations=[angle, -angle - π], scatter_kwargs...)

# vertical lines
x = collect(0.38:0.01:0.75)
yup = @. -(x - 1) * (x + 1) + 0.1 * (1 / (x - 0.3) - 1)
ydown = @. -(x - 1) * (x + 1) - 0.1 * (1 / (x - 0.3) - 1)
lines!(ax_graphical[:maj], x, yup; line_kwargs...)
lines!(ax_graphical[:maj], -x, yup; line_kwargs...)
lines!(ax_graphical[:maj], x, yup; line_kwargs...)
lines!(ax_graphical[:maj], -x, yup; line_kwargs...)
lines!(ax_graphical[:maj], x[ydown .> 0.05], ydown[ydown .> 0.05]; line_kwargs...)
lines!(ax_graphical[:maj], -x[ydown .> 0.05], ydown[ydown .> 0.05]; line_kwargs...)

ind = 8
angle = atan(yup[ind] - yup[ind - 1], x[ind] - x[ind - 1])
scatter!(ax_graphical[:maj], [x[ind], -x[ind]], [yup[ind], yup[ind]];
         rotations=[angle, -angle - π], scatter_kwargs...)

ind = 10
angle = atan(ydown[ind] - ydown[ind - 1], x[ind] - x[ind - 1])
scatter!(ax_graphical[:maj], [x[ind], -x[ind]], [ydown[ind], ydown[ind]];
         rotations=[angle, -angle - π], scatter_kwargs...)

# stable points
scatter!(ax_graphical[:maj], [0], [1]; marker=:circle,
         color=:white, markersize=10)
scatter!(ax_graphical[:maj], [0], [1]; marker='◑',
         color=:black,
         markersize=12)
scatter!(ax_graphical[:maj], [-0.4, 0, 0.4], [0, 0, 0]; marker=:circle,
         color=:white, markersize=10)
scatter!(ax_graphical[:maj], [-0.4, 0, 0.4], [0, 0, 0]; marker='○',
         color=:black, markersize=10)

# PROPORTIONAL VOTING
# vertical lines
lines!(ax_graphical[:prop], [-0.5, -0.5], [0.1, 2]; line_kwargs...)
lines!(ax_graphical[:prop], [0.5, 0.5], [0.1, 2]; line_kwargs...)

scatter!(ax_graphical[:prop], [-0.5, 0.5, -0.5, 0.5],
         [1.5, 1.5, 0.4, 0.4];
         rotations=[-π / 2, -π / 2, π / 2, π / 2], scatter_kwargs...)

x = [-0.5, 0, 0.5]
scatter!(ax_graphical[:prop], x, -(x .- 1) .* (x .+ 1); marker=:circle,
         color=:black, markersize=10)
scatter!(ax_graphical[:prop], x, [0, 0, 0]; marker=:circle,
         color=:white, markersize=10)
scatter!(ax_graphical[:prop], x, [0, 0, 0]; marker=:circle,
         color=:white, markersize=10)
scatter!(ax_graphical[:prop], x, [0, 0, 0]; marker='○',
         color=:black, markersize=10)

leg = Legend(fig[1, 2], axis[:prop_rtr]; merge=true, orientation=:vertical, labelsize=12)

colgap!(ax_gridpos, 15)
rowgap!(ax_gridpos, 15)

display(fig)

if prompt
    prompt_for_save(filename, fig)
end