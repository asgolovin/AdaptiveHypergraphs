using AdaptiveHypergraphs
const AH = AdaptiveHypergraphs

include("figure_plotting_tools.jl")

# ======================================================
# -------------------- INPUT ---------------------------

input_folder = joinpath(projectdir(),
                        "data/thesis/run_2023-02-23_13-32-26_motifs_D4")
data_folder = DataFolder(input_folder, :simulation)

num_nodes = 10000
max_size = 4
triples = [OrderTwoMotif((3, 0), (1, 0), (1, 0)),
           OrderTwoMotif((1, 1), (1, 0), (0, 1)),
           OrderTwoMotif((0, 3), (1, 0), (1, 0)),
           OrderTwoMotif((0, 2), (1, 0), (0, 2))]
mc_motifs = filter(x -> AH.order(x) <= 1, all_motifs(max_size))

triple_labels = Dict(OrderTwoMotif((0, 2), (1, 0), (0, 2)) => "[B²|A|B²]",
                     OrderTwoMotif((1, 1), (1, 0), (0, 1)) => "[AB|A|B]",
                     OrderTwoMotif((0, 3), (1, 0), (1, 0)) => "[B³|A|A]",
                     OrderTwoMotif((3, 0), (1, 0), (1, 0)) => "[A³|A|A]")
mc_labels = Dict(OrderTwoMotif((0, 2), (1, 0), (0, 2)) => "1/2 [AB²][AB²]/[A]",
                 OrderTwoMotif((1, 1), (1, 0), (0, 1)) => "2 [A²B][AB]/[A]",
                 OrderTwoMotif((0, 3), (1, 0), (1, 0)) => "2 [AB³][A²]/[A]",
                 OrderTwoMotif((3, 0), (1, 0), (1, 0)) => "8 [A⁴][A²]/[A]")

node_colormap = :RdYlGn_6
hyperedge_colormap = :thermal
misc_colormap = :Egypt

num_triples = length(triples)
linecolors = get(tum_scheme, 1:4, (1, 4))

skip_points = 20

prompt = true

filename = "./figures/thesis/moment_closure.pdf"
ylabel = "# of motifs"
ylims = (0, 5000)
xticks = 0:20:100
yticks = 0:1000:5000

# ------------------------------------------------------
# ======================================================

# create figure
size_inches = (5, 3)
size_pt = 1.5 * 72 .* size_inches
fig = Figure(; resolution=size_pt)

ax = Axis(fig[1, 1];
          titlesize=15,
          title="Moment closure approximation",
          xticklabelsize=12,
          yticklabelsize=12,
          xlabelsize=13,
          ylabelsize=13,
          titlefont="Latin Modern Roman")
ax.xlabel = L"time $t$"
ax.ylabel = ylabel
ax.xticks = xticks
ax.yticks = yticks
xlims!(ax, (0, 100))
ylims!(ax, ylims)

for (batch_folder, batch_num) in data_folder
    for (run_folder, run_num) in batch_folder
        if run_num != 1
            continue
        end

        mc_motif_values = Dict(motif => [] for motif in mc_motifs)
        for motif in mc_motifs
            if AH.order(motif) == 0
                if motif.state == AH.A
                    indices, values = load_data("state_count_A.csv", run_folder)
                    indices = indices[1:skip_points:end]
                    values = values[1:skip_points:end]
                else
                    indices, values = load_data("state_count_B.csv", run_folder)
                    indices = indices[1:skip_points:end]
                    values = values[1:skip_points:end]
                end
            else
                indices, values = load_data("motif_count_$motif.csv", run_folder)
                indices = indices[1:skip_points:end]
                values = values[1:skip_points:end]
            end
            mc_motif_values[motif] = values
        end

        for (i, triple_motif) in enumerate(triples)
            indices, triple = load_data("motif_count_$triple_motif.csv", run_folder)
            indices = indices[1:skip_points:end]
            triple = triple[1:skip_points:end]

            left_motif = triple_motif.left_motif
            right_motif = triple_motif.right_motif
            int_state = triple_motif.int.A > 0 ? AH.A : AH.B
            intersection = OrderZeroMotif(int_state)

            # compute the combinatorical prefactor
            if int_state == AH.A
                left_count = left_motif.A
                right_count = right_motif.A
            else
                left_count = left_motif.B
                right_count = right_motif.B
            end
            issymmetrical = left_motif == right_motif ? 0.5 : 1
            prefactor = issymmetrical * left_count * right_count

            closure = prefactor .* mc_motif_values[left_motif] .*
                      mc_motif_values[right_motif] ./
                      mc_motif_values[intersection]

            println("$triple_motif: prefactor: $prefactor")

            lines!(ax, indices, triple; label=triple_labels[triple_motif],
                   color=linecolors[i], linewidth=1.5)
            lines!(ax, indices, closure;
                   label=mc_labels[triple_motif],
                   color=linecolors[i] * 0.6, linewidth=1.5)
        end
    end
end

group_lines = [LineElement(; color=linecolors[i], linewidth=2) for i in 1:4]
group_lines_dark = [LineElement(; color=linecolors[i] * 0.6, linewidth=2) for i in 1:4]
group_labels = [triple_labels[triples[i]] for i in 1:4]
group_labels_dark = [mc_labels[triples[i]] for i in 1:4]

leg = Legend(fig[1, 2],
             collect(collect.(zip(group_lines, group_lines_dark))),
             collect(collect.(zip(group_labels, group_labels_dark))),
             ["", "", "", ""];
             orientation=:vertical, labelsize=14)
leg.titleposition = :left
leg.gridshalign = :left
leg.groupgap = 13
display(fig)

if prompt
    prompt_for_save(filename, fig; pt_per_unit=0.666)
end