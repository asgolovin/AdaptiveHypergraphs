using AdaptiveHypergraphs
const AH = AdaptiveHypergraphs

const Label = AH.Label

include("figure_plotting_tools.jl")

# ======================================================
# -------------------- INPUT ---------------------------

input_folder = joinpath(projectdir(),
                        "final_data/run_2023-02-23_13-32-26_motifs_D4")
data_folder = DataFolder(input_folder, :simulation)

num_nodes = 10000
max_size = 4
triples = [Label("[A3|A|A]"), Label("[AB|A|B]"), Label("[B3|A|A]"), Label("[B2|A|B2]")]
mc_labels = filter(x -> AH.order(x) <= 1, all_labels(max_size))

triple_latex_labels = Dict(Label("[B2|A|B2]") => "[B²|A|B²]",
                           Label("[AB|A|B]") => "[AB|A|B]",
                           Label("[B3|A|A]") => "[B³|A|A]",
                           Label("[A3|A|A]") => "[A³|A|A]")
mc_latex_labels = Dict(Label("[B2|A|B2]") => "1/2 [AB²][AB²]/[A]",
                       Label("[AB|A|B]") => "2 [A²B][AB]/[A]",
                       Label("[B3|A|A]") => "2 [AB³][A²]/[A]",
                       Label("[A3|A|A]") => "8 [A⁴][A²]/[A]")

node_colormap = :RdYlGn_6
hyperedge_colormap = :thermal
misc_colormap = :Egypt

num_triples = length(triples)
linecolors = get(tum_scheme, 1:4, (1, 4))

skip_points = 20

prompt = true

filename = "./figures/moment_closure.pdf"
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

        mc_motifs = Dict(label => [] for label in mc_labels)
        for label in mc_labels
            if AH.order(label) == 0
                if label.left[AH.A] == 1
                    indices, values = load_data("state_count_A.csv", run_folder)
                    indices = indices[1:skip_points:end]
                    values = values[1:skip_points:end]
                else
                    indices, values = load_data("state_count_B.csv", run_folder)
                    indices = indices[1:skip_points:end]
                    values = values[1:skip_points:end]
                end
            else
                indices, values = load_data("motif_count_$label.csv", run_folder)
                indices = indices[1:skip_points:end]
                values = values[1:skip_points:end]
            end
            mc_motifs[label] = values
        end

        for (i, triple_label) in enumerate(triples)
            indices, triple = load_data("motif_count_$triple_label.csv", run_folder)
            indices = indices[1:skip_points:end]
            triple = triple[1:skip_points:end]

            left_label = Label(triple_label.left_total)
            right_label = Label(triple_label.right_total)
            intersection = Label(triple_label.int)

            # compute the combinatorical prefactor
            int_state = triple_label.int_state
            left_count = triple_label.left_total[int_state]
            right_count = triple_label.right_total[int_state]
            issymmetrical = left_label == right_label ? 0.5 : 1
            prefactor = issymmetrical * left_count * right_count

            closure = prefactor .* mc_motifs[left_label] .* mc_motifs[right_label] ./
                      mc_motifs[intersection]

            println("$triple_label: prefactor: $prefactor")

            lines!(ax, indices, triple; label=triple_latex_labels[triple_label],
                   color=linecolors[i], linewidth=1.5)
            lines!(ax, indices, closure;
                   label=mc_latex_labels[triple_label],
                   color=linecolors[i] * 0.6, linewidth=1.5)
        end
    end
end

group_lines = [LineElement(; color=linecolors[i], linewidth=2) for i in 1:4]
group_lines_dark = [LineElement(; color=linecolors[i] * 0.6, linewidth=2) for i in 1:4]
group_labels = [triple_latex_labels[triples[i]] for i in 1:4]
group_labels_dark = [mc_latex_labels[triples[i]] for i in 1:4]

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