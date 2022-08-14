using DataStructures
using Random
using Distributions

export AbstractModel, DiscrModel, ContinuousModel, step!


abstract type AbstractModel{P <: PropagationRule, A <: AdaptivityRule} end


struct DiscrModel{P <: PropagationRule, A <: AdaptivityRule} <: AbstractModel{P, A}
    network::HyperNetwork
    propagation_rule::P
    adaptivity_rule::A
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


mutable struct ContinuousModel{P <: PropagationRule, A <: AdaptivityRule} <: AbstractModel{P, A}
    network::HyperNetwork
    event_queue::PriorityQueue
    # for every hyperedge, stores the time of the next event
    next_event_time::Dict{Int64, Float64}
    propagation_rule::P
    adaptivity_rule::A
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

    propagation_distr = Exponential(propagation_rate)
    adaptivity_distr = Exponential(adaptivity_rate)

    # add all active hyperedges to the queue
    next_event_time = Dict{Int64, Float64}()
    for hyperedge in get_hyperedges(network)
        next_event_time[hyperedge] = Inf
        if is_active(network, hyperedge)
            time = rand(propagation_distr)
            event = Event(hyperedge, time, propagate)
            next_event_time[hyperedge] = 0.
            enqueue!(event_queue, event, event.time)
        end
    end

    current_time = 0.

    ContinuousModel{P, A}(network,
                    event_queue,
                    next_event_time,
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

    # Don't do anything if there are no events anymore
    if length(model.event_queue) == 0
        return false
    end

    event = dequeue!(model.event_queue)
    model.current_time = event.time

    source_hyperedge = event.hyperedge

    # Should not be necessary... theoretically
    if !is_active(network, source_hyperedge)
        println("A non-active hyperedge was selected by an event. Something went wrong!")
        return false
    end

    # hyperedges which share a common node with the source hyperedge
    neighboring_hyperedges = Dict{Int64, Bool}()
    # nodes which changed their state after the event
    affected_nodes = Dict{Int64, NamedTuple{(:before, :after), Tuple{State, State}}}()

    if event.action == propagate
        # record whether the neighbors are active or not
        for node in get_nodes(network, source_hyperedge)
            for h in get_hyperedges(network, node)
                if !(h == source_hyperedge || h in keys(neighboring_hyperedges))
                    neighboring_hyperedges[h] = is_active(network, h)
                end
            end
        end

        println("Executing propagation rule")
        affected_nodes = propagate!(network, propagation_rule, source_hyperedge)
        
    elseif event.action == adapt
        # since an adaptivity rule can target any hyperedge, we need to record the whole state of the network.
        # This can be written in a faster way if this ever becomes a bottleneck
        for h in get_hyperedges(network)
            if h != source_hyperedge
                neighboring_hyperedges[h] = is_active(network, h)
            end
        end

        println("Executing adaptivity rule")
        affected_nodes = adapt!(network, adaptivity_rule, source_hyperedge)
    end

    # Add events for the source hyperedge
    # if the source hyperedeg is still active, the Poisson process is restarted
    if is_active(network, source_hyperedge)
        _add_event!(model, source_hyperedge, propagate)
    else # if the source hyperedge was switched off, all future events are removed
        _remove_hyperedge_events!(model.event_queue, source_hyperedge)
    end

    # Add events for the neighboring hyperedges
    for neighbor in keys(neighboring_hyperedges)
        # if any node of the neighbor was affected
        # TODO: doesn't work for adaptation!
        if any(keys(affected_nodes) .∈ Ref(get_nodes(network, neighbor)))

            active_before = neighboring_hyperedges[neighbor]
            active_after = is_active(network, neighbor)

            # on -> off
            if active_before == true && active_after == false
                _remove_hyperedge_events!(model.event_queue, neighbor)
            
            # off -> on
            elseif active_before == false && active_after == true
                _add_event!(model, neighbor, propagate)

            end # in all other cases, nothing happens
        end
    end
    
    return length(affected_nodes) > 0
end


"""
    _remove_hyperedge_events!(queue::PriorityQueue, hyperedge::Integer)

Remove all events belonging to the hyperedge `hyperedge`
"""
function _remove_hyperedge_events!(queue::PriorityQueue, hyperedge::Integer)
    for (event, _) in queue
        if event.hyperedge == hyperedge
            delete!(queue, event)
        end
    end
end

function _add_event!(model::ContinuousModel, hyperedge::Integer, event_type::EventTypes)
    distr = event_type == propagate ? model.propagation_distr : model.adaptivity_distr
    event_time = model.current_time + rand(distr)
    event = Event(hyperedge, event_time, event_type)
    enqueue!(model.event_queue, event, event.time)
end