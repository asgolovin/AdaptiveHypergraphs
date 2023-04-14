using AdaptiveHypergraphs
using MPI

default_input = "../input/default_network.jl"

# command-line interface
if length(ARGS) > 0
    input_file = ARGS[1]
    include(input_file)
else
    input_file = default_input
end

include(input_file)

start_simulation(params)