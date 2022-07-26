using DataStructures

include("../Model.jl")

struct ModelObservable{P <: PropagationRule, A <: AdaptivityRule}
    model::Observable{TS{P, A} <: AbstractModel}
    state_hist::Dict{State, CircularBuffer{Int64}}
end

function step!(mo::ModelObservable)
    model, state_hist = mo.model, mo.state_hist
    step!(model)
end