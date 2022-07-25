using SimpleHypergraphs
using GraphPlot
using ColorSchemes
using GLMakie

# ====================================================================================
# ---------------------------------- TEST CODE ---------------------------------------

include("simulation/Model.jl")
include("presentation/HypergraphPlot.jl")

n = 100
network = Observable(HyperNetwork(n, 0.5))
build_RSC_hg!(network[], (2n, n ÷ 3, 1, 1))

majority_rule = MajorityRule()
rewiring_rule = RewiringRule(0.5)

model = DiscrModel{MajorityRule, RewiringRule}(network[],
                                                            majority_rule,
                                                            rewiring_rule)

f = Figure()
display(f)

hypergraphplot(f[1, 1], network)

for i = 1:1000
    network_changed = step!(model)
    if network_changed
        network[] = network[]
    end
    sleep(0.1)
end