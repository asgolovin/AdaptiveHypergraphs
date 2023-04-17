include("figure_plotting_tools.jl")

# ======================================================
# -------------------- INPUT ---------------------------

input_folder_rtr = joinpath(projectdir(),
                            "data/thesis/run_2023-02-17_18-00-07_simple_graph_rtr")
input_folder_rts = joinpath(projectdir(),
                            "data/thesis/run_2023-02-17_18-00-08_simple_graph_rts")

num_nodes = 10000
max_size = 4
hyperedge_colormap = :thermal
linecolors = get(colorschemes[hyperedge_colormap], [1, 2.5, 3], (1, max_size))
skip_initial = 200
skip_mid = 30
skip_last = 30000
prop_rule = :rts

if prop_rule == :rtr
    input_folder = input_folder_rtr
    filename = "./figures/thesis/slow_manifold_graph_rtr.pdf"
else
    input_folder = input_folder_rts
    filename = "./figures/thesis/slow_manifold_graph_rts.pdf"
end

# ------------------------------------------------------
# ======================================================

# create figure
size_inches = (5, 3)
size_pt = 1.5 .* 72 .* size_inches
fig = Figure(; resolution=size_pt)
ax = Axis(fig[1, 1];
          titlesize=15,
          xlabel=L"magnetization $m$",
          ylabel=L"density of active links $\rho_{AB}$",
          title=L"Slow manifold plot $$",
          xticklabelsize=12,
          yticklabelsize=12,
          xlabelsize=13,
          ylabelsize=13)

magnetization_data = []
density_data = []

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

        filename = "state_count_A.csv"
        path = joinpath(rundir, filename)
        df = CSV.read(path, DataFrame; delim=", ")

        state_count_A = df.value
        magnetization = 2 .* (state_count_A ./ num_nodes) .- 1
        if batch_num == 1
            magnetization = vcat(magnetization[1:skip_initial],
                                 magnetization[skip_initial:skip_mid:(end - skip_last)])
            println(length(magnetization))
        end

        filename = "active_hyperedge_count_2.csv"
        path = joinpath(rundir, filename)
        df = CSV.read(path, DataFrame; delim=", ")

        active_hyperedeges = df.value
        if batch_num == 1
            active_hyperedeges = vcat(active_hyperedeges[1:skip_initial],
                                      active_hyperedeges[skip_initial:skip_mid:(end - skip_last)])
        end
        active_density = active_hyperedeges ./ num_nodes

        if batch_num == 1
            append!(magnetization_data, magnetization[skip_initial:end])
            append!(density_data, active_density[skip_initial:end])
        end

        linecolor = batch_num == 1 ? linecolors[1] : linecolors[3]
        label = batch_num == 1 ? L"$m_0 = 0$" : L"$m_0 = -0.8$"

        # plot one trajectory in a darker color and with a label
        if run_num == 4
            lines!(ax, magnetization, active_density; color=linecolor,
                   label=label, linewidth=0.7)
        else
            lines!(ax, magnetization, active_density; color=(linecolor, 0.3),
                   linewidth=0.7)
        end
    end
end

xlims!(ax, (-1, 1))
ylims!(ax, (0, 1))

p, cov = fit_parabola(magnetization_data, density_data, false)
x = collect(-1:0.01:1)
lines!(ax, x, p.(x); color=:red, label="polynomial fit")

leg = Legend(fig[1, 2], ax; merge=true, orientation=:vertical, labelsize=12)

display(fig)

ans = Base.prompt("Press space and then twice Enter")
options = ["no", "yes"]
menu = RadioMenu(options)

choice = request("Save the figure to file $filename?", menu)
if choice == 2
    save(filename, fig; pt_per_unit=0.666)
end