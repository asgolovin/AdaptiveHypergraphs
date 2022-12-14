using Random

export AdaptivityRule, RewireToRandom, RewireToSame, adapt!

"""
    AdaptivityRule

Determines how the edges of the graph change with time. 
"""
abstract type AdaptivityRule end

"""
    RewireToRandom <: AdaptivityRule

A randomly chosen node from the hyperedge leaves the hyperedge and connects to 
either a different randomly chosen node or a different hyperedge. 
"""
struct RewireToRandom <: AdaptivityRule end

"""
    ConflictAvoiding <: AdaptivityRule

Similar to the RewireToRandom rule, but the selected node only connects to nodes or hyperedges in the same state.
"""
struct RewireToSame <: AdaptivityRule end

"""
A randomly chosen node from the hyperedge leaves the hyperedge and connects to 
either a different randomly chosen node or a different hyperedge. 

Return a list of modified hyperedges.
"""
function adapt!(network::HyperNetwork, adaptivity_rule::RewireToRandom, hyperedge::Integer)
    selected_node = rand(get_nodes(network, hyperedge))

    function hyperedge_conditions(h)
        size = get_hyperedge_size(network, h)
        max_size = get_max_size(network)
        nodes = get_nodes(network, h)
        return h != hyperedge && # has to be a different hyperedge
               size < max_size && # the size should not exceed the limit
               selected_node ∉ nodes # no degenerate hyperdeges
    end

    function node_conditions(node)
        return node != selected_node
    end

    hyperedge_candidates, node_candidates = _rejection_sampling(network,
                                                                hyperedge_conditions,
                                                                node_conditions)

    affected_hyperedges = _rewire_to_candidate!(network,
                                                hyperedge,
                                                selected_node,
                                                hyperedge_candidates,
                                                node_candidates)

    return affected_hyperedges
end

"""
Similar to the RewireToRandom rule, but the selected node only connects to nodes or hyperedges in the same state.

Return a list of modified hyperedges.
"""
function adapt!(network::HyperNetwork, adaptivity_rule::RewireToSame,
                hyperedge::Integer)
    selected_node = rand(get_nodes(network, hyperedge))
    required_state = get_state(network, selected_node)

    function hyperedge_conditions(h)
        size = get_hyperedge_size(network, h)
        max_size = get_max_size(network)
        nodes = get_nodes(network, h)
        return !is_active(network, h) &&  # all nodes of the hyperedge have to be in the same state
               get_state(network, nodes[1]) == required_state && # the state has to be equal to the state of the node
               size < max_size && # the size should not exceed the limit
               selected_node ∉ nodes # no degenerate hyperdeges
    end

    function node_conditions(node)
        return node != selected_node &&
               get_state(network, node) == required_state
    end

    hyperedge_candidates, node_candidates = _rejection_sampling(network,
                                                                hyperedge_conditions,
                                                                node_conditions)

    affected_hyperedges = _rewire_to_candidate!(network,
                                                hyperedge,
                                                selected_node,
                                                hyperedge_candidates,
                                                node_candidates)
    return affected_hyperedges
end

function _rejection_sampling(network, hyperedge_conditions, node_conditions)
    candidate_found = false
    hyperedge_candidates = []
    node_candidates = []
    while !candidate_found
        r = rand(1:(get_num_hyperedges(network) + get_num_nodes(network)))
        if r <= get_num_hyperedges(network)
            new_hyperedge = get_hyperedges(network)[r]
            if hyperedge_conditions(new_hyperedge)
                candidate_found = true
                hyperedge_candidates = [new_hyperedge]
                node_candidates = []
            end
        else
            node = r - get_num_hyperedges(network)
            if node_conditions(node)
                candidate_found = true
                hyperedge_candidates = []
                node_candidates = [node]
            end
        end
    end
    return hyperedge_candidates, node_candidates
end

function _rewire_to_candidate!(network, old_hyperedge, selected_node, hyperedge_candidates,
                               node_candidates)
    num_candidates = length(hyperedge_candidates) + length(node_candidates)

    if num_candidates == 0
        println("The network has no other nodes or hyperedges that node $selected_node can connect with.")
        return []
    end

    candidate_id = rand(1:num_candidates)

    #print("Node $selected_node was disconnected from nodes $(get_nodes(network, old_hyperedge)) ")

    if candidate_id <= length(hyperedge_candidates)
        new_hyperedge = hyperedge_candidates[candidate_id]
        # print("and connected to nodes $(get_nodes(network, new_hyperedge))\n")
        remove_node!(network, selected_node, old_hyperedge)
        include_node!(network, selected_node, new_hyperedge)
    else
        new_node = node_candidates[candidate_id - length(hyperedge_candidates)]
        # print("and connected to node $new_node\n")
        remove_node!(network, selected_node, old_hyperedge)
        new_hyperedge = add_hyperedge!(network, [selected_node, new_node])
    end

    affected_hyperedges = [new_hyperedge]
    return affected_hyperedges
end