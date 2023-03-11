# first results before improvements: (commit 90b401a68ce) -> with a sysimage
# time from start of main.jl script to empty figure: 1:50 min -> 0:34 min
# time to first print: 0:36 min, 2:26 total -> 0:30 min, 1:04 total
# output of `@time AdaptiveHypergraphs`: 
# >> 79.965480 seconds (140.74 M allocations: 8.807 GiB, 2.12% gc time, 5.40% compilation time: 66% of which was recompilation)

# include("data/Motif.jl")  -> 0.49 sec
# include("data/Network.jl") -> 0.11 sec

# include("simulation/AdaptivityRule.jl")   -> 0.008 sec
# include("simulation/PropagationRule.jl")  -> 0.004 sec
# include("simulation/Model.jl")            -> 0.66 sec
# include("simulation/InputParams.jl")    -> 2.00 sec
# include("simulation/MomentExpansion.jl")  -> 9.05 sec!

# include("presentation/HypergraphPlot.jl") -> 22.65 sec
# include("presentation/Measurements.jl")   -> 0.05 sec
# include("presentation/ModelObservable.jl") -> 0.52 sec
# include("presentation/Panels.jl")         -> 0.11 sec
# include("presentation/Dashboard.jl")      -> 0.12 sec
# include("presentation/SimulationBatch.jl") -> 1.00 sec

# execute @time for the first time to precompile it
@time 2 + 2

tic = time()

println("GLMakie:")
@time using GLMakie # 0.011 sec

println("DifferentialEquations:")
@time using DifferentialEquations # 0.00022 sec

println("SimpleHypergraphs:")
@time using SimpleHypergraphs # 0.00022 sec

println("Label.jl:")
@time include("../src/data/Motif.jl") # 0.87 sec

println("Network.jl:")
@time include("../src/data/Network.jl") # 0.14 sec

println("AdaptivityRule.jl:")
@time include("../src/simulation/AdaptivityRule.jl") # 0.007

println("PropagationRule.jl:")
@time include("../src/simulation/PropagationRule.jl") # 0.005

println("Model.jl:")
@time include("../src/simulation/Model.jl") # 0.085 sec

println("InputParams.jl:")
@time include("../src/simulation/InputParams.jl") # 1.97 !

println("MomentExpansion.jl:")
@time include("../src/simulation/MomentExpansion.jl") # 0.084 sec

println("HypergraphPlot.jl:")
@time include("../src/presentation/HypergraphPlot.jl") # 1.01

println("Measurements.jl:")
@time include("../src/presentation/Measurements.jl") # 0.19

println("ModelObservable.jl:")
@time include("../src/presentation/ModelObservable.jl") # 0.92

println("Panels.jl:")
@time include("../src/presentation/Panels.jl") # 0.21

println("Dashboard.jl:")
@time include("../src/presentation/Dashboard.jl") # 0.29

println("SimulationBatch.jl:")
@time include("../src/presentation/SimulationBatch.jl") # 1.01 -> 0.029

toc = time()

println("Total time: $(toc - tic)")