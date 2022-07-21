module AdaptiveHypergraphs

# Use HyperNetX to plot the graph
using PyCall
using Conda
# Conda.runconda(`install matplotlib --yes`)
# Conda.runconda(`install networkx --yes`)
# run(`$(PyCall.python) -m pip install 'decorator>=5.0.9'`)
# run(`$(PyCall.python) -m pip install hypernetx`)

using SimpleHypergraphs
using PyPlot
using GraphPlot
using ColorSchemes

# use python matplotlib backend
pygui(true)

# ====================================================================================
# ---------------------------------- TEST CODE ---------------------------------------

include("simulation/TimeStepper.jl")

n = 4
network = HyperNetwork(n, 0.5)
build_RSC_hg!(network, (4, 3))

majority_rule = MajorityRule()
rewiring_rule = RewiringRule(0.7)

time_stepper = DiscrTimeStepper{MajorityRule, RewiringRule}(network, majority_rule, rewiring_rule)

for i = 1:10
    step!(time_stepper)
    println(network.hg)

    #SimpleHypergraphs.draw(network.hg, HyperNetX; width=5, height=5, no_border=true)
    #show()
    #readline()
end



end