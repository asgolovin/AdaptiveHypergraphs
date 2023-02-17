using AdaptiveHypergraphs
using MPI

default_input = "../input/large_network.jl"

# command-line interface
if length(ARGS) > 0
    input_file = ARGS[1]
    try
        include(input_file)
    catch SystemError
        throw(ArgumentError("The file $input_file doesn't exist. Please enter a valid input file."))
    end
else
    input_file = default_input
end

include(input_file)

start_simulation(params)