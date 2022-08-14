using DataStructures
using Random
using Distributions

export AbstractModel, DiscrModel, ContinuousModel, step!


abstract type AbstractModel{P <: PropagationRule, A <: AdaptivityRule} end


struct DiscrModel{P <: PropagationRule, A <: AdaptivityRule} <: AbstractModel{P, A}
    network::HyperNetwork
    propagation_rule::PropagationRule
    adaptivity_rule::AdaptivityRule
    propagation_prob::Real
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
    hyperedge = rand(get_hyperedges(network))

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


struct ContinuousModel{P <: PropagationRule, A <: AdaptivityRule} <: AbstractModel{P, A}
    network::HyperNetwork
    event_queue::PriorityQueue
    propagation_rule::PropagationRule
    adaptivity_rule::AdaptivityRule
    propagation_distr::Distribution
    adaptivity_distr::Distribution
    current_time::Real
end


function ContinuousModel{P, A}(network::HyperNetwork,
                               propagation_rule::P,
                               adaptivity_rule::A, 
                               propagation_rate::Real,
                               adaptivity_rate::Real) where {P <: PropagationRule, A <: AdaptivityRule}
    event_queue = PriorityQueue()
    he_meta = (:change_of_status_time => Inf)

    propagation_distr = Exponential(propagation_rate)
    adaptivity_distr = Exponential(adaptivity_rate)

    # add all active hyperedges to the queue
    for hyperedge in get_hyperedges(network)
        set_hyperedge_meta!(network, hyperedge, he_meta)
        if is_active(network, hyperedge)
            time = rand(propagation_distr)
            event = Event(hyperedge, time, propagate)
            set_hyperedge_meta!(network, hyperedge, :change_of_status_time, 0.)
            enqueue!(event_queue, event, event.time)
        end
    end

    current_time = 0.

    ContinuousModel(network,
                    event_queue,
                    propagation_rule,
                    adaptivity_rule,
                    propagation_distr,
                    adaptivity_distr,
                    current_time)
end

@enum EventTypes propagate=0 adapt=1

struct Event
    hyperedge::Integer
    time::Real
    action::EventTypes
end


function step!(model::ContinuousModel)
    network = model.network
    propagation_rule = model.propagation_rule
    adaptivity_rule = model.adaptivity_rule

    if length(model.event_queue) == 0
        return nothing
    end

    event = dequeue!(model.event_queue)
    model.current_time = event.time

    source_hyperedge = event.hyperdege

    if !is_active(network, source_hyperedge)
        return nothing
    end

    if event.action == propagate
        # record whether the neighbors are active or not
        neighboring_hyperedges = Dict{Int64, Bool}()
        for node in get_nodes(network, source_hyperedge)
            for h in get_hyperedges(network, node)
                if !(h == source_hyperedge || h in keys(neighboring_hyperedges))
                    neighboring_hyperedges[h] = is_active(network, h)
                end
            end
        end

        println("Executing propagation rule")
        affected_nodes = propagate!(network, propagation_rule, source_hyperedge)

        # check if the source hyperedge was modified 
        # if the source hyperedge was modified but is still active or was not modified at all
        if is_active(network, source_hyperedge)
            # add a new event for this hyperedge
            event_time = model.current_time + rand(model.propagation_distr)
            event = Event(source_hyperedge, event_time, propagate)
            enqueue!(model.event_queue, event)

        else # if the source hyperedge was switched off
            # remove all future events for this hyperedge
            for event in model.event_queue
                if event.hyperdege == source_hyperedge
                    delete!(model.event_queue, event)
                end
            end
        end

        # check if any neighboring hyperedges were modified
        for neighbor in keys(neighboring_hyperedges)
            # if any node of the neighbor was affected
            if any(keys(affected_nodes) .∈ Ref(get_nodes(network, neighbor)))
                active_before = neighboring_hyperedges[neighbor]
                active_after = is_active(network, neighbor)
                if active_before == true && active_after == false
                    # remove all future events for this hyperdege
                    for event in model.event_queue
                        if event.hyperdege == neighbor
                            delete!(model.event_queue, event)
                        end
                    end

                elseif active_before == false && active_after == true
                    # add new event for this hyperedge
                    event_time = model.current_time + rand(model.propagation_distr)
                    event = Event(neighbor, event_time, propagate)
                    enqueue!(model.event_queue, event)

                end # in all other cases, nothing happens
            end
        end

    elseif event.action == adapt
        # TODO
        # println("Executing adaptivity rule")
        # adapt!(network, adaptivity_rule, source_hyperedge)

        # process changed hyperedges

    end
end