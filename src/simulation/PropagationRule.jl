using StatsBase

export PropagationRule, MajorityRule, propagate!

"""
    PropagationRule

Determines how the states of the graph change with time. 
"""
abstract type PropagationRule end


"""
All nodes within the hyperedge take the opinion of the majority. 

In case of a tie, a decision is made randomly. 
"""
struct MajorityRule <: PropagationRule end


function propagate!(network::HyperNetwork, majority_rule::MajorityRule, hyperedge::Integer)
    nodes = get_nodes(network, hyperedge)
    node_to_state = get_node_to_state_dict(network, hyperedge)
    state_count = countmap(values(node_to_state))
    # number of votes for the majority opinion
    max_count = maximum(values(state_count))

    # track which nodes were affected by the change
    affected_nodes = Int64[]

    # get all opinions with the same maximum number of votes
    max_states = [pair.first for pair in state_count if pair.second == max_count]
    if length(max_states) == 1
        majority_state = max_states[1]
    else
        majority_state = rand(max_states)
    end
    
    for node in nodes
        prev_state = get_state(network, node)
        if prev_state != majority_state
            push!(affected_nodes, node)
            set_state!(network, node, majority_state)
        end
    end

    println("Nodes $(nodes) in hyperedge $hyperedge were set to $majority_state")

    return affected_nodes
end