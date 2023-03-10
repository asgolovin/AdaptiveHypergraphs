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
function adapt!(network::HyperNetwork, adaptivity_rule::RewireToRandom, hyperedge::Int64)
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
                hyperedge::Int64)
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

    candidate_id, candidate_type = _rejection_sampling(network,
                                                       hyperedge_conditions,
                                                       node_conditions)

    affected_hyperedges = _rewire_to_candidate!(network,
                                                hyperedge,
                                                selected_node,
                                                candidate_id,
                                                candidate_type)
    return affected_hyperedges
end

"""
    _rejection_sampling(network, hyperedge_conditions, node_conditions)

Select candidates for rewiring until a candidate which fulfills all conditions is found.

# Arguments
- `network`: the hypergraph
- `hyperedge_conditions`: a function which takes the index of a hyperedge and returns true if the hyperedge can be rewired to, false otherwise
- `node_conditions`: same, but for the node index

# Return
- `candidade_id::Int64`: the id of the node or hyperedge. Return `nothing` if the sampling was not successful. 
- `candidate_type::Symbol`: `:node` or `:hyperedge` depending on which was chosen. Return `nothing` if the sampling was not successful. 
"""
function _rejection_sampling(network::HyperNetwork,
                             hyperedge_conditions::Function,
                             node_conditions::Function)
    num_samples = 0
    max_samples = 100 # stop sampling if the max number of tries is exceeded
    while num_samples < max_samples
        num_samples += 1
        r = rand(1:(get_num_hyperedges(network) + get_num_nodes(network)))
        if r <= get_num_hyperedges(network)
            new_hyperedge = get_hyperedges(network)[r]
            if hyperedge_conditions(new_hyperedge)
                candidate_id = new_hyperedge
                candidite_type = :hyperedge
                return candidate_id, candidite_type
            end
        else
            node = r - get_num_hyperedges(network)
            if node_conditions(node)
                candidate_id = node
                candidite_type = :node
                return candidate_id, candidite_type
            end
        end
    end
    return nothing, nothing
end

function _rewire_to_candidate!(network::HyperNetwork, old_hyperedge::Int64,
                               selected_node::Int64, candidate_id::Int64,
                               candidate_type::Symbol)
    if candidate_type == :hyperedge
        new_hyperedge = candidate_id
        remove_node!(network, selected_node, old_hyperedge)
        include_node!(network, selected_node, new_hyperedge)
    elseif candidate_type == :node
        new_node = candidate_id
        remove_node!(network, selected_node, old_hyperedge)
        new_hyperedge = add_hyperedge!(network, [selected_node, new_node])
    else
        throw(ArgumentError("candidate_type has to be either :hyperedge or :node, not $candidate_type"))
    end

    affected_hyperedges = [new_hyperedge]
    return affected_hyperedges
end