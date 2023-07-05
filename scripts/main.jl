using AdaptiveHypergraphs
using MPI

input_file = "../input/default_network.jl"

# command-line interface
if length(ARGS) > 0
    input_file = ARGS[1]
end

include(input_file)

start_simulation(params)