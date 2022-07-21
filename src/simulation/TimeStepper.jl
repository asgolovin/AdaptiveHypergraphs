using DataStructures
using Random

include("../data/Network.jl")
include("AdaptivityRule.jl")
include("PropagationRule.jl")


abstract type AbstractTimeStepper{P <: PropagationRule, A <: AdaptivityRule} end


struct DiscrTimeStepper{P <: PropagationRule, A <: AdaptivityRule} <: AbstractTimeStepper{P, A}
    network::HyperNetwork
    propagation_rule::PropagationRule
    adaptivity_rule::AdaptivityRule
    #...
end


struct ContinuousTimeStepper{P <: PropagationRule, A <: AdaptivityRule} <: AbstractTimeStepper{P, A}
    network::HyperNetwork
    event_queue::PriorityQueue
    propagation_rule::PropagationRule
    adaptivity_rule::AdaptivityRule
    #...
end


"""
Advances the dynamics of the network by one step. 
"""
function step!(time_stepper::DiscrTimeStepper)
    network = time_stepper.network
    propagation_rule = time_stepper.propagation_rule
    adaptivity_rule = time_stepper.adaptivity_rule

    # choose a random hyperedge
    hyperedge = rand(1:get_num_hyperedges(network))

    # do nothing if the hyperedge connects vertices with the same state
    if !is_active(network, hyperedge)
        println("Doing nothing")
        return nothing
    end

    p = rand()
    if p < adaptivity_rule.rewiring_prob
        println("Executing adaptivity rule")
        adapt!(network, adaptivity_rule, hyperedge)
    else
        println("Executing propagation rule")
        propagate!(network, propagation_rule, hyperedge)
    end

    return nothing
end

function step!(time_stepper::ContinuousTimeStepper)

    # return τ
end