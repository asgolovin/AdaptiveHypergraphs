using SimpleHypergraphs
using GraphPlot
using ColorSchemes
using GLMakie

# ====================================================================================
# ---------------------------------- TEST CODE ---------------------------------------

include("simulation/TimeStepper.jl")
include("presentation/HypergraphPlot.jl")

n = 30
network = Observable(HyperNetwork(n, 0.5))
build_RSC_hg!(network[], (2n, n ÷ 3, 1, 1))

majority_rule = MajorityRule()
rewiring_rule = RewiringRule(0.5)

time_stepper = DiscrTimeStepper{MajorityRule, RewiringRule}(network[],
                                                            majority_rule,
                                                            rewiring_rule)

f = Figure()
display(f)

hypergraphplot(f[1, 1], network)

while true
    step!(time_stepper)
    network[] = network[]
    sleep(0.1)
end