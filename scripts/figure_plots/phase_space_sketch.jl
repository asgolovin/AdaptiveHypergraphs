interactive = true

if interactive
    using GLMakie
    GLMakie.activate!()
else
    include("figure_plotting_tools.jl")

    filename = "figures/phase_plot_sketch.pdf"
    prompt = false
end

function draw_phase_plot!(ax, p, parabola)
    xs = collect(-1.0:0.01:1.0)
    lines_x = [-0.5, 0, 0.5]

    ymax = parabola(0, 0)

    parabola_values = @lift parabola(xs, $p)
    down_triangle_positions = @lift ymax * 2 .- 0.5 .* (ymax * 2 .- parabola(lines_x, $p))

    up_triangle_positions = lift(p) do p
        if p < p_c
            return 0.5 .* parabola(lines_x, p)
        else
            return fill(NaN, size(lines_x))
        end
    end
    parabola_circle_positions = @lift parabola(lines_x, $p)

    lines!(ax, xs, zero(xs); line_kwargs...)
    lines!(ax, xs, parabola_values; line_kwargs...)

    for x in lines_x
        lines!(ax, [x, x], [ymax * 1.9, ymax * 0.1]; line_kwargs...)
    end
    scatter!(ax, lines_x, down_triangle_positions; rotations=-π / 2, triangle_kwargs...)
    scatter!(ax, lines_x, up_triangle_positions; rotations=π / 2, triangle_kwargs...)

    scatter!(ax, lines_x, zero(lines_x); empty_circle_background_kwargs...)
    scatter!(ax, lines_x, zero(lines_x); empty_circle_kwargs...)
    scatter!(ax, lines_x, parabola_circle_positions; circle_kwargs...)
    scatter!(ax, [-1, 1], [0, 0]; circle_kwargs...)

    ylims!(ax, (-ymax * 0.1, ymax * 2))
    return ax
end

if interactive
    fig = Figure()
    sg = SliderGrid(fig[2, 1],
                    (label=L"p", range=0:0.01:1, startvalue=0.3))
    p = sg.sliders[1].value

    ax = Axis(fig[1, 1]; titlesize=24, xlabelsize=22,
              ylabelsize=22)

    ax.title = "phase space sketch"

    hideydecorations!(ax; label=false)
    ax.xlabel = "Magnetization"
    ax.ylabel = "Density of active edges"
    ax.xticks = [-1, 0, 1]
else
    fig = create_figure(:small, 1.7)
    ax_dict = Dict(:low => Axis(fig[1, 1]; titlesize=14),
                   :high => Axis(fig[1, 2]; titlesize=14))
    for ax in values(ax_dict)
        hidexdecorations!(ax; label=false, ticks=false, ticklabels=false)
        hideydecorations!(ax; label=false)
        ax.xlabel = L"m"
        ax.xticks = [-1, 0, 1]
    end
    ax_dict[:low].ylabel = L"\rho"

    ax_dict[:low].title = L"p\,<\,p_c"
    ax_dict[:high].title = L"p\,>\,p_c"

    text!(ax_dict[:low], 0.05, 0.98; text="A", font=Makie.theme(nothing, :fonts).bold[],
          space=:relative, align=(:left, :top),
          textsize=16)
    text!(ax_dict[:high], 0.05, 0.98; text="B", font=Makie.theme(nothing, :fonts).bold[],
          space=:relative, align=(:left, :top),
          textsize=16)

    colgap!(fig.layout, 15)
end

# drawing settings
if interactive
    line_kwargs = Dict(:color => :black, :linewidth => 2)
    triangle_kwargs = Dict(:color => :black, :marker => :rtriangle, :markersize => 20)
    circle_kwargs = Dict(:color => :black, :marker => :circle, :markersize => 30)
    empty_circle_background_kwargs = Dict(:color => :white, :marker => :circle,
                                          :markersize => 30)
    empty_circle_kwargs = Dict(:color => :black, :marker => '○', :markersize => 30)
else
    line_kwargs = Dict(:color => :black, :linewidth => 1)
    triangle_kwargs = Dict(:color => :black, :marker => :rtriangle)
    circle_kwargs = Dict(:color => :black, :marker => :circle)
    empty_circle_background_kwargs = Dict(:color => :white, :marker => :circle)
    empty_circle_kwargs = Dict(:color => :black, :marker => '○')
end

const p_c = 0.7

if interactive
    parabola(x, p) = max.(-x .^ 2 .+ 1.0 .- 1.5 .* p, 0.0)
    draw_phase_plot!(ax, p, parabola)
else
    low_p = Observable(0)
    high_p = Observable(1)
    parabola(x, p) = max.(-x .^ 2 .+ 1.0 .- 1.5 .* p, 0.0)

    # fit the parabola to the data
    data_folder = joinpath(projectdir(),
                           "data/run_2023-05-31_13-53-20_simple_graph_phase_space_sketch")
    run_folder = Dict(:low => joinpath(data_folder, "batch_001", "run_001"),
                      :high => joinpath(data_folder, "batch_002", "run_001"))

    magnetization = Dict(:low => Float64, :high => Float64[])
    active_edge_density = Dict(:low => Float64, :high => Float64[])

    for panel in [:low, :high]
        _, state_count = load_data(joinpath(run_folder[panel], "state_count_A.csv"),
                                   run_folder[panel])
        num_nodes = 1000
        magnetization[panel] = (state_count ./ num_nodes) .* 2 .- 1

        _, active_edges = load_data(joinpath(run_folder[panel],
                                             "active_hyperedge_count_2.csv"),
                                    run_folder[panel])
        active_edge_density[panel] = active_edges ./ num_nodes

        lines!(ax_dict[panel], magnetization[panel], active_edge_density[panel];
               linewidth=0.7,
               color=hyperedge_linecolors(4)[3])
    end

    fitted_parabola, _ = fit_parabola(magnetization[:low], active_edge_density[:low])
    parabola(x, p) = p < p_c ? fitted_parabola.(x) : zero(x)

    draw_phase_plot!(ax_dict[:low], low_p, parabola)
    draw_phase_plot!(ax_dict[:high], high_p, parabola)
end

display(fig)

if !interactive && prompt
    prompt_for_save(filename, fig)
end