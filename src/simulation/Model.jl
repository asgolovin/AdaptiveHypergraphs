using DataStructures
using Random

export AbstractModel, DiscrModel, ContinuousModel, step!


abstract type AbstractModel{P <: PropagationRule, A <: AdaptivityRule} end


struct DiscrModel{P <: PropagationRule, A <: AdaptivityRule} <: AbstractModel{P, A}
    network::HyperNetwork
    propagation_rule::PropagationRule
    adaptivity_rule::AdaptivityRule
    propagation_prob::Real
end


struct ContinuousModel{P <: PropagationRule, A <: AdaptivityRule} <: AbstractModel{P, A}
    network::HyperNetwork
    event_queue::PriorityQueue
    propagation_rule::PropagationRule
    adaptivity_rule::AdaptivityRule
    propagation_rate::Real
    adaptivity_rate::Real
end


"""
Advances the dynamics of the network by one step. 

Returns true if the network has changed, false otherwise. 
"""
function step!(model::DiscrModel)
    network = model.network
    propagation_rule = model.propagation_rule
    adaptivity_rule = model.adaptivity_rule

    # choose a random hyperedge
    hyperedge = rand(1:get_num_hyperedges(network))

    # do nothing if the hyperedge connects vertices with the same state
    if !is_active(network, hyperedge)
        return false
    end

    p = rand()
    if p < model.propagation_prob
        println("Executing propagation rule")
        propagate!(network, propagation_rule, hyperedge)
    else
        println("Executing adaptivity rule")
        adapt!(network, adaptivity_rule, hyperedge)
    end

    return true
end

# TODO
function step!(model::ContinuousModel)
    network = model.network
    propagation_rule = model.propagation_rule
    adaptivity_rule = model.adaptivity_rule

    
end