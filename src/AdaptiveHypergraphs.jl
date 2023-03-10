module AdaptiveHypergraphs

WITH_DISPLAY = get(ENV, "DISPLAY", "") != ""

include("data/Label.jl")
include("data/Network.jl")

include("simulation/AdaptivityRule.jl")
include("simulation/PropagationRule.jl")
include("simulation/Model.jl")
include("simulation/InputParams.jl")
include("simulation/MomentExpansion.jl")

@static if WITH_DISPLAY
    include("presentation/HypergraphPlot.jl")
end
include("presentation/Measurements.jl")
include("presentation/ModelObservable.jl")
@static if WITH_DISPLAY
    include("presentation/Panels.jl")
end
include("presentation/NinjaDashboard.jl")
@static if WITH_DISPLAY
    include("presentation/Dashboard.jl")
else
    struct Dashboard <: AbstractDashboard end
    function Dashboard(::Any, ::Any; save_folder::Any)
        return nothing
    end
end
include("presentation/SimulationBatch.jl")

end