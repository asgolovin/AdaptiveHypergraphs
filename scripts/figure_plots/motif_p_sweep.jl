using StatsBase

using AdaptiveHypergraphs
const AH = AdaptiveHypergraphs

include("figure_plotting_tools.jl")

# ======================================================
# -------------------- INPUT ---------------------------

draw_error_bars = false

#"results/run_2024-01-18_12-07-59_motifs_p_sweep_maj_rts") # long simulation
folder_maj_rts = joinpath(projectdir(),
                          "results/run_2024-01-21_15-42-22_motifs_p_sweep_maj_rts") # short simulation 
folder_maj_rtr = joinpath(projectdir(),
                          "results/run_2024-01-21_15-36-25_motifs_p_sweep_maj_rtr")
folder_prop_rts = joinpath(projectdir(),
                           "results/run_2024-01-21_15-50-56_motifs_p_sweep_prop_rts")
folder_prop_rtr = joinpath(projectdir(),
                           "results/run_2024-01-21_15-53-41_motifs_p_sweep_prop_rtr")

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
linecolors = hyperedge_linecolors(max_size)

prompt = false

filename = "./figures/motifs_p_sweep.pdf"

xlabel = L"$p$"
ylabel = L"[\Xi]/[\Xi]_\text{closure}"
xticks = 0:0.2:1
yticks = 0:1:6
xlims = (-0.02, 1.02)
ylims = (-0.02, 6)

max_Δp = 1
max_Δm = 0.06

ABA_triples = [OrderTwoMotif((1, 0), (0, 1), (1, 0)), # [A(B)A]
               OrderTwoMotif((2, 0), (0, 1), (1, 0)), # [A2(B)A]
               OrderTwoMotif((2, 0), (0, 1), (2, 0)), # [A2(B)A2]
               OrderTwoMotif((3, 0), (0, 1), (1, 0))] # [A3(B)A]

BBA_triples = [OrderTwoMotif((0, 1), (0, 1), (1, 0)), # [B(B)A]
               OrderTwoMotif((0, 2), (0, 1), (1, 0)), # [B2(B)A]
               OrderTwoMotif((0, 2), (0, 1), (2, 0))] # [B2(B)A2]

mixed_triples = [OrderTwoMotif((1, 1), (0, 1), (1, 0)), # [AB(B)A]
                 OrderTwoMotif((1, 1), (0, 1), (1, 1))] # [AB(B)AB]

#triples = vcat(ABA_triples, BBA_triples, mixed_triples)

triples = [motif for motif in all_motifs(max_size) if AH.order(motif) == 2]

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

        mu = motif.int.A
        nu = motif.int.B

        if (mu + nu) > 1
            continue
        end

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

axis = add_four_rules_axes!(ax_gridpos, titles,
                            xlabel, ylabel,
                            xticks, yticks,
                            xlims, ylims)
#axis.title = "Ratio of simulated triples to the moment closure"

for panel in panels
    hlines!(axis[panel], [1.0]; color=:gray, linestyle=:dash, linewidth=1)

    p_values = Float64[]
    ratios = Dict(motif => Float64[] for motif in triples)
    errors = Dict(motif => Float64[] for motif in triples)

    data_folder = DataFolder(input_folders[panel], :simulation)

    for (batch_folder, batch_num) in data_folder
        params = load_params(joinpath(batch_folder.folder, "input_params.json"))
        p = params.model_params.adaptivity_prob
        push!(p_values, p)

        for motif in triples
            ratio, error = motif_ratio(batch_folder, motif)
            push!(ratios[motif], ratio)
            push!(errors[motif], error)
        end
    end

    indices = sortperm(p_values)
    p_values = p_values[indices]
    for motif in triples
        ratios[motif] = ratios[motif][indices]
        errors[motif] = errors[motif][indices]
    end

    for motif in triples
        #scatter!(axis[panel], p_values, ratios[motif];
        #        markersize=10, marker=:xcross, label="$motif")
        lines!(axis[panel], p_values, ratios[motif]; linewidth=1.5, label="$motif")
        if draw_error_bars
            errorbars!(axis[panel], p_values, ratios[motif], errors[motif];
                       whiskerwidth=7, linewidth=0.5)
        end
    end
end

# for motif in BBA_triples
#     scatter!(axis, p_values, ratios[motif];
#              markersize=10, marker=:rect, label="$motif")
#     if draw_error_bars
#         errorbars!(axis, p_values, ratios[motif], errors[motif];
#                    whiskerwidth=7, linewidth=0.5)
#     end
# end

# for motif in mixed_triples
#     scatter!(axis, p_values, ratios[motif];
#              markersize=10, marker=:circle, label="$motif")
#     if draw_error_bars
#         errorbars!(axis, p_values, ratios[motif], errors[motif];
#                    whiskerwidth=7, linewidth=0.5)
#     end
# end

Legend(fig[1, 2], axis[:maj_rtr]; merge=true, orientation=:vertical, labelsize=12)

display(fig)

if prompt
    prompt_for_save(filename, fig)
end
