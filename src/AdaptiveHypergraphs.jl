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

dashboard = Dashboard(model, 100; plot_hypergraph=true, is_interactive=false)
run!(dashboard, 1000, 10)
# record!(dashboard, "test_record", 100, 10)
end