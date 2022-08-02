module AdaptiveHypergraphs

using SimpleHypergraphs
using ColorSchemes
using GLMakie

include("data/Network.jl")

include("simulation/AdaptivityRule.jl")
include("simulation/PropagationRule.jl")
include("simulation/Model.jl")

include("presentation/HypergraphPlot.jl")
include("presentation/ModelObservable.jl")
include("presentation/Dashboard.jl")

# ====================================================================================
# ---------------------------------- TEST CODE ---------------------------------------

n = 30
network = Observable(HyperNetwork(n, 0.5))
build_RSC_hg!(network[], (2n, n ÷ 3, 1, 1))

majority_rule = MajorityRule()
rewiring_rule = RewiringRule(0.5)

model = DiscrModel{MajorityRule, RewiringRule}(network[],
                                               majority_rule,
                                               rewiring_rule)

create_dashboard(model, 20; plot_hypergraph=true, interactivity=false)

# record(fig, "results\\first_test.mp4", 1:100, framerate = 3, compression = 1) do i
#     network_changed = step!(model)
#     if network_changed
#         network[] = network[]
#     end
#     ax.title = "t = $i"
# end

# for i = 1:10
#     network_changed = step!(model)
#     if network_changed
#         network[] = network[]
#     end
#     sleep(0.1)
# end

end