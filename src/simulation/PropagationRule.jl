using StatsBase

export PropagationRule, MajorityVoting, ProportionalVoting, propagate!

"""
    PropagationRule

Determines how the states of the graph change with time. 
"""
abstract type PropagationRule end

"""
    MajorityVoting <: PropagationRule

All nodes within the hyperedge become convinced of the opinion of the majority. 
In case of a tie, a decision is made randomly. 
"""
struct MajorityVoting <: PropagationRule end

"""
ProportionalVoting <: PropagationRule

An opinion is chosen with a probability that is proportional to 
the share of nodes with this opinion. All nodes within the hyperedge
become convinced of this opinion. 
"""
struct ProportionalVoting <: PropagationRule end

"""
    propagate!(network::HyperNetwork, majority_voting::MajorityVoting,
               hyperedge::Integer)

All nodes within the hyperedge become convinced of the opinion of the majority. 
In case of a tie, a decision is made randomly. 

Returns a list of modified nodes.
"""
function propagate!(network::HyperNetwork, majority_voting::MajorityVoting,
                    hyperedge::Integer)
    nodes = get_nodes(network, hyperedge)
    node_to_state = get_node_to_state_dict(network, hyperedge)
    state_count = countmap(values(node_to_state))
    # number of votes for the majority opinion
    max_count = maximum(values(state_count))

    # get all opinions with the same maximum number of votes
    max_states = [pair.first for pair in state_count if pair.second == max_count]
    if length(max_states) == 1
        majority_state = max_states[1]
    else
        majority_state = rand(max_states)
    end

    affected_nodes = _convince!(network, hyperedge, majority_state)

    return affected_nodes
end

"""
    propagate!(network::HyperNetwork, proportional_voting::ProportionalVoting,
               hyperedge::Integer)

An opinion is chosen with a probability that is proportional to 
the share of nodes with this opinion. All nodes within the hyperedge
become convinced of this opinion. 

Returns a list of modified nodes.
"""
function propagate!(network::HyperNetwork, proportional_voting::ProportionalVoting,
                    hyperedge::Integer)
    nodes = get_nodes(network, hyperedge)
    trendsetter = rand(nodes)

    # choosing one opinion proportional to the share of "votes" is equivalent 
    # to choosing one random node and taking his opinion
    chosen_state = get_state(network, trendsetter)
    affected_nodes = _convince!(network, hyperedge, chosen_state)

    return affected_nodes
end

function _convince!(network::HyperNetwork, hyperedge::Integer, state::State)
    affected_nodes = Int64[]
    nodes = get_nodes(network, hyperedge)
    for node in nodes
        prev_state = get_state(network, node)
        if prev_state != state
            push!(affected_nodes, node)
            set_state!(network, node, state)
        end
    end
    println("Nodes $affected_nodes in hyperedge $hyperedge were set to $state")
    return affected_nodes
end