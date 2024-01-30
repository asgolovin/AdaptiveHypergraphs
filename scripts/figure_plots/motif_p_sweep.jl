using StatsBase
using ColorSchemes

using AdaptiveHypergraphs
const AH = AdaptiveHypergraphs

include("figure_plotting_tools.jl")

# ======================================================
# -------------------- INPUT ---------------------------

draw_error_bars = false

folder_maj_rts = joinpath(projectdir(),
                          "results/run_2024-01-28_13-31-36_motifs_p_sweep_maj_rts")
folder_maj_rtr = joinpath(projectdir(),
                          "results/run_2024-01-28_17-19-32_motifs_p_sweep_maj_rtr")
folder_prop_rts = joinpath(projectdir(),
                           "results/run_2024-01-28_12-03-05_motifs_p_sweep_prop_rts")
folder_prop_rtr = joinpath(projectdir(),
                           "results/run_2024-01-28_10-40-31_motifs_p_sweep_prop_rtr")

panels = [:maj_rts, :maj_rtr, :prop_rts, :prop_rtr]
titles = Dict(:maj_rts => "Majority voting,\nrewire-to-same",
              :maj_rtr => "Majority voting,\nrewire-to-random",
              :prop_rts => "Proportional voting\nrewire-to-same",
              :prop_rtr => "Proportional voting\nrewire-to-random")

input_folders = Dict(:maj_rts => folder_maj_rts,
                     :maj_rtr => folder_maj_rtr,
                     :prop_rts => folder_prop_rts,
                     :prop_rtr => folder_prop_rtr)

num_nodes = 1000
max_size = 4

prompt = false

filename = "./figures/motifs_p_sweep.pdf"

xlabel = L"$p$"
ylabel = L"[\Xi]/[\Xi]_\text{closure}"
xticks = 0:0.2:1
rts_yticks = 0:10:100
rtr_yticks = 0:0.5:3
xlims = (-0.02, 1.02)
rts_ylims = (0.0, 55)
rtr_ylims = (0.0, 3)

max_Δp = 1
max_Δm = 0.06

groups = [:ABA, :BBA, :mixed]

triples = Dict(:ABA => [OrderTwoMotif((1, 0), (0, 1), (1, 0)), # [A(B)A]
                        OrderTwoMotif((2, 0), (0, 1), (1, 0)), # [A2(B)A]
                        OrderTwoMotif((2, 0), (0, 1), (2, 0)), # [A2(B)A2]
                        OrderTwoMotif((3, 0), (0, 1), (1, 0))], # [A3(B)A]
               :BBA => [OrderTwoMotif((0, 1), (0, 1), (1, 0)), # [B(B)A]
                        OrderTwoMotif((0, 2), (0, 1), (1, 0)), # [B2(B)A]
                        OrderTwoMotif((0, 2), (0, 1), (2, 0))], # [B2(B)A2]
               :mixed => [OrderTwoMotif((1, 1), (0, 1), (1, 0)), # [AB(B)A]
                          OrderTwoMotif((1, 1), (0, 1), (0, 1)), # [AB(B)B]
                          OrderTwoMotif((1, 1), (0, 1), (1, 1)), # [AB(B)AB]
                          OrderTwoMotif((2, 1), (0, 1), (1, 1))]) # [A2B(B)AB]

num_triples = Dict(group => length(triples[group]) for group in groups)

labels = Dict(group => [pretty_label(motif) for motif in triples[group]]
              for group in groups)

colormaps = Dict(:ABA => :Oranges, :BBA => :Greens, :mixed => :Purples)

colors = Dict(group => get(colorschemes[colormaps[group]], collect(num_triples[group]:-1:1),
                           (-1, num_triples[group]))
              for group in groups)

# ------------------------------------------------------
# ======================================================

function motif_ratio(batch_folder, motif)
    ratio_vector = Float64[]

    for (run_folder, _) in batch_folder
        simulated_motif_vector = load_data("motif_count_$(motif).csv", run_folder)[2]

        if maximum(simulated_motif_vector) < 10
            continue
        end

        left_motif = motif.left_motif
        right_motif = motif.right_motif

        A_count = load_data("state_count_A.csv", run_folder)[2]
        B_count = load_data("state_count_B.csv", run_folder)[2]

        left_motif_count = load_data("motif_count_$(left_motif).csv", run_folder)[2]
        right_motif_count = load_data("motif_count_$(right_motif).csv", run_folder)[2]

        # ignore datapoints where one of the counts is zero
        mask = simulated_motif_vector .!= 0

        simulated_motif_vector = simulated_motif_vector[mask]
        left_motif_count = left_motif_count[mask]
        right_motif_count = right_motif_count[mask]
        A_count = A_count[mask]
        B_count = B_count[mask]

        mu = motif.int.A
        nu = motif.int.B

        mc_motif_vector = left_motif_count .* right_motif_count ./ binomial.(A_count, mu) ./
                          binomial.(B_count, nu)
        mc_motif_vector *= moment_closure_coeff(motif)

        push!(ratio_vector, mean(simulated_motif_vector ./ mc_motif_vector))
    end

    if length(ratio_vector) < 3
        return NaN, NaN
    end

    ratio = mean(ratio_vector)
    error = std(ratio_vector)

    return ratio, error
end

# create figure
fig = create_figure(:large, 3 / 2)
ax_gridpos = fig[1, 1] = GridLayout()

panels = [:maj_rts, :maj_rtr, :prop_rts, :prop_rtr]
axis = Dict()
for (i, panel) in enumerate(panels)
    row = (i + 1) ÷ 2
    col = mod1(i, 2)
    axis[panel] = Axis(ax_gridpos[row, col])
    axis[panel].title = titles[panel]

    if row == 1
        hidexdecorations!(axis[panel]; grid=false)
    else
        axis[panel].xlabel = xlabel
    end

    axis[panel].xticks = xticks

    if panel == :maj_rts || panel == :prop_rts
        axis[panel].yticks = rts_yticks
        ylims!(axis[panel], rts_ylims)
        axis[panel].ylabel = ylabel
    else
        axis[panel].yticks = rtr_yticks
        ylims!(axis[panel], rtr_ylims)
    end
    xlims!(axis[panel], xlims)
end
colgap!(ax_gridpos, 15)
rowgap!(ax_gridpos, 15)

for panel in panels
    hlines!(axis[panel], [1.0]; color=:gray, linestyle=:dash, linewidth=1)

    p_values = Float64[]
    ratios = Dict(group => Dict(motif => Float64[] for motif in triples[group])
                  for group in groups)
    errors = Dict(group => Dict(motif => Float64[] for motif in triples[group])
                  for group in groups)

    data_folder = DataFolder(input_folders[panel], :simulation)

    for (batch_folder, batch_num) in data_folder
        params = load_params(joinpath(batch_folder.folder, "input_params.json"))
        p = params.model_params.adaptivity_prob
        push!(p_values, p)

        for group in groups
            for motif in triples[group]
                ratio, error = motif_ratio(batch_folder, motif)
                push!(ratios[group][motif], ratio)
                push!(errors[group][motif], error)
            end
        end
    end

    indices = sortperm(p_values)
    p_values = p_values[indices]
    for group in groups
        for motif in triples[group]
            ratios[group][motif] = ratios[group][motif][indices]
            errors[group][motif] = errors[group][motif][indices]
        end
    end

    for group in groups
        for (i, motif) in enumerate(triples[group])
            scatter!(axis[panel], p_values, ratios[group][motif];
                     markersize=10, marker=:xcross, color=colors[group][i])
            lines!(axis[panel], p_values, ratios[group][motif]; linewidth=1,
                   color=colors[group][i])
            if draw_error_bars
                errorbars!(axis[panel], p_values, ratios[group][motif],
                           errors[group][motif];
                           whiskerwidth=7, linewidth=0.5, color=colors[group][i])
            end
        end
    end
end

# draw the legend
elements = Dict(group => [[LineElement(; color=color, linewidth=1),
                           MarkerElement(; color=color, marker='x', markersize=10)]
                          for color in colors[group]]
                for group in groups)

leg = Legend(fig[1, 2],
             [elements[:ABA], elements[:BBA], elements[:mixed]],
             [labels[:ABA], labels[:BBA], labels[:mixed]],
             ["Aᵃ ( B ) Aᵐ motifs:", "Bᵇ ( B ) Aᵐ motifs:", "Mixed motifs:"])

leg.gridshalign = :left
leg.gridsvalign = :top

display(fig)

if prompt
    prompt_for_save(filename, fig)
end
