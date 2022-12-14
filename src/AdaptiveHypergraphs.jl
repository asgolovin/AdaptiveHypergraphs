module AdaptiveHypergraphs

include("data/Label.jl")
include("data/Network.jl")

include("simulation/AdaptivityRule.jl")
include("simulation/PropagationRule.jl")
include("simulation/Model.jl")

include("presentation/InputParams.jl")
include("presentation/HypergraphPlot.jl")
include("presentation/Measurements.jl")
include("presentation/ModelObservable.jl")
include("presentation/Panels.jl")
include("presentation/Dashboard.jl")
include("presentation/SimulationBatch.jl")

end