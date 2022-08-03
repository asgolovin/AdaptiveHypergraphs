using DataStructures
using Random

export AbstractModel, DiscrModel, ContinuousModel, step!


abstract type AbstractModel{P <: PropagationRule, A <: AdaptivityRule} end


struct DiscrModel{P <: PropagationRule, A <: AdaptivityRule} <: AbstractModel{P, A}
    network::HyperNetwork
    propagation_rule::PropagationRule
    adaptivity_rule::AdaptivityRule
    #...
end


struct ContinuousModel{P <: PropagationRule, A <: AdaptivityRule} <: AbstractModel{P, A}
    network::HyperNetwork
    event_queue::PriorityQueue
    propagation_rule::PropagationRule
    adaptivity_rule::AdaptivityRule
    #...
end


"""
Advances the dynamics of the network by one step. 

Returns true if the netwrok has changed, false otherwise. 
"""
function step!(model::DiscrModel)
    network = model.network
    propagation_rule = model.propagation_rule
    adaptivity_rule = model.adaptivity_rule

    # choose a random hyperedge
    hyperedge = rand(1:get_num_hyperedges(network))

    # do nothing if the hyperedge connects vertices with the same state
    if !is_active(network, hyperedge)
        println("Doing nothing")
        return false
    end

    p = rand()
    if p < adaptivity_rule.rewiring_prob
        println("Executing adaptivity rule")
        adapt!(network, adaptivity_rule, hyperedge)
    else
        println("Executing propagation rule")
        propagate!(network, propagation_rule, hyperedge)
    end

    return true
end

# TODO
function step!(model::ContinuousModel)

    # return τ
end