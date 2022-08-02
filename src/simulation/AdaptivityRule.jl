using Random

"""
    AdaptivityRule

Determines how the edges of the graph change with time. 
"""
abstract type AdaptivityRule end

struct RewiringRule <: AdaptivityRule
    rewiring_prob::Real
end


"""
A randomly chosen node from the hyperedge leaves the hyperedge and connects to a different node with the same state. 
Roughly like in "An adaptive voter model on simplical complexes"
"""
function adapt!(network::HyperNetwork, adaptivity_rule::RewiringRule, hyperedge::Integer)
    nodes = get_nodes(network, hyperedge)
    selected_node = rand(nodes)
    required_state = get_state(network, selected_node)

    # find a new random node with the same state 
    state_dict = get_node_to_state_dict(network)
    filter!(pair -> (pair.second == required_state) && (pair.first != selected_node), state_dict)

    if length(state_dict) == 0
        println("The network has no other nodes in state $required_state that node $selected_node can connect with.")
        return nothing
    end
    
    new_node = rand(keys(state_dict))

    remove_node_from_hyperedge!(network, selected_node, hyperedge)
    add_hyperedge!(network, [selected_node, new_node])

    println("Node $selected_node was disconnected from hyperedge $hyperedge and connected with node $new_node")
    return nothing
end