using Random

export AdaptivityRule, RewiringRule, ConflictAvoiding, adapt!

"""
    AdaptivityRule

Determines how the edges of the graph change with time. 
"""
abstract type AdaptivityRule end

"""
    RewiringRule <: AdaptivityRule

A randomly chosen node from the hyperedge leaves the hyperedge and connects to 
either a different randomly chosen node or a different hyperedge. 
"""
struct RewiringRule <: AdaptivityRule end

"""
    ConflictAvoiding <: AdaptivityRule

Similar to the RewiringRule rule, but the selected node only connects to nodes or hyperedges in the same state.
"""
struct ConflictAvoiding <: AdaptivityRule end

"""
A randomly chosen node from the hyperedge leaves the hyperedge and connects to 
either a different randomly chosen node or a different hyperedge. 

Returns a list of modified hyperedges.
"""
function adapt!(network::HyperNetwork, adaptivity_rule::RewiringRule, hyperedge::Integer)
    selected_node = rand(get_nodes(network, hyperedge))

    # find hyperedges in the same state 
    hyperedge_candidates = get_hyperedges(network)
    function conditions(h)
        size = get_hyperedge_size(network, h)
        max_size = get_max_hyperedge_size(network)
        return h != hyperedge && # has to be a different hyperedge
               size < max_size # the size should not exceed the limit
    end
    filter!(conditions, hyperedge_candidates)

    # find nodes in the same state
    node_candidates = get_nodes(network)
    filter!(node -> node != selected_node, node_candidates)

    affected_hyperedges = _rewire_to_candidate!(network,
                                                hyperedge,
                                                selected_node,
                                                hyperedge_candidates,
                                                node_candidates)

    return affected_hyperedges
end

"""
Similar to the RewiringRule rule, but the selected node only connects to nodes or hyperedges in the same state.

Returns a list of modified hyperedges.
"""
function adapt!(network::HyperNetwork, adaptivity_rule::ConflictAvoiding,
                hyperedge::Integer)
    selected_node = rand(get_nodes(network, hyperedge))
    state_dict = get_node_to_state_dict(network)
    required_state = state_dict[selected_node]

    # find hyperedges in the same state 
    hyperedge_candidates = get_hyperedges(network)
    function conditions(h)
        size = get_hyperedge_size(network, h)
        max_size = get_max_hyperedge_size(network)
        nodes = get_nodes(network, h)
        return !is_active(network, h) &&  # all nodes of the hyperedge have to be in the same state
               state_dict[nodes[1]] == required_state && # the state has to be equal to the state of the node
               size < max_size # the size should not exceed the limit
    end
    filter!(conditions, hyperedge_candidates)

    # find nodes in the same state
    node_candidates = get_nodes(network)
    filter!(node -> (state_dict[node] == required_state) && (node != selected_node),
            node_candidates)

    affected_hyperedges = _rewire_to_candidate!(network,
                                                hyperedge,
                                                selected_node,
                                                hyperedge_candidates,
                                                node_candidates)
    return affected_hyperedges
end

function _rewire_to_candidate!(network, old_hyperedge, selected_node, hyperedge_candidates,
                               node_candidates)
    num_candidates = length(hyperedge_candidates) + length(node_candidates)

    if num_candidates == 0
        println("The network has no other nodes or hyperedges that node $selected_node can connect with.")
        return []
    end

    candidate_id = rand(1:num_candidates)

    print("Node $selected_node was disconnected from nodes $(get_nodes(network, old_hyperedge)) ")

    if candidate_id <= length(hyperedge_candidates)
        new_hyperedge = hyperedge_candidates[candidate_id]
        print("and connected to nodes $(get_nodes(network, new_hyperedge))\n")
        remove_node_from_hyperedge!(network, selected_node, old_hyperedge)
        add_node_to_hyperedge!(network, selected_node, new_hyperedge)
    else
        new_node = node_candidates[candidate_id - length(hyperedge_candidates)]
        print("and connected to node $new_node\n")
        remove_node_from_hyperedge!(network, selected_node, old_hyperedge)
        new_hyperedge = add_hyperedge!(network, [selected_node, new_node])
    end

    affected_hyperedges = [new_hyperedge]
    return affected_hyperedges
end