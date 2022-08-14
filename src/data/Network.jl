using SimpleHypergraphs
using StatsBase
using Random
using Graphs
using Combinatorics

export State, HyperNetwork, build_regular_hg!, build_RSC_hg!

@enum State::Bool I=false S=true

"""
    HyperNetwork

A higher-order network where every node can be in a particular state (SIS or SIR).

The topology of the network is represented by a hypergraph: a graph, where every hyperedge can connect multiple vertices, and not just two. 

Other than the hypergraph and the states, the struct keeps track of any statistics related to the system such as the number of infected nodes. 
"""
mutable struct HyperNetwork
    # The underlying hypergraph
    hg::Hypergraph{Bool, State}
    # Number of nodes in a particular state
    state_dist::Dict{State, Integer}
    # Number of hyperedges of a particular size
    hyperedge_dist::Dict{Int64, Int64}
    # A bijective map where the vector indices correspond to the indices of the columns 
    # in the incidence matrix of the hypergraph (matrix IDs or MIDs) and the values 
    # to unique IDs (UIDs). The UIDs start with 1. 
    # Only the UIDs are exposed outside of the class; the MIDs should 
    # only be used to communicate with SimpleHypergraphs. 
    # This mapping is important because the native hyperedge indices from 
    # SimpleHypergraphs are not preserved when a hyperedge is deleted. 
    hyperedge_uid::Vector{Int64}
    # current highest hyperedge UID
    max_hyperedge_uid::Integer
end


# ====================================================================================
# ------------------------------- CONSTRUCTORS ---------------------------------------


"""
    HyperNetwork(n::Integer, node_state::Vector{Union{Nothing, State}})

Create an empty network with n nodes and no hyperedges.
`node_state` denotes the state of each node. 
"""
function HyperNetwork(n::Integer,
                      node_state::Vector{Union{Nothing, State}})
    @assert length(node_state) == n
    matrix = Matrix{Union{Nothing, Bool}}(nothing, (n, 0))
    hg = Hypergraph{Bool, State}(matrix; v_meta=node_state)
    state_dist = countmap(node_state)
    hyperedge_dist = Dict(2 => 0)
    hyperedge_uid = Vector{Int64}()
    HyperNetwork(hg, state_dist, hyperedge_dist, hyperedge_uid, 0)
end


"""
    HyperNetwork(n::Integer)

Create an empty network with n nodes and no hyperedges.
All nodes are suseptible. 
"""
function HyperNetwork(n::Integer)
    node_state = Vector{Union{Nothing, State}}(nothing, n)
    fill!(node_state, S)
    matrix = Matrix{Union{Nothing, Bool}}(nothing, (n, 0))
    hg = Hypergraph{Bool, State}(matrix; v_meta=node_state)
    hyperedge_dist = Dict(2 => 0)
    hyperedge_uid = Vector{Int64}()
    HyperNetwork(hg, Dict(S => n, I => 0), hyperedge_dist, hyperedge_uid, 0)
end


"""
    HyperNetwork(n::Integer, p0::AbstractFloat)

Create an empty network with n nodes and no hyperedges.
Each node is infected with probability p0. 
"""
function HyperNetwork(n::Integer, p0::AbstractFloat)
    node_state = Vector{Union{Nothing, State}}(nothing, n)
    for i in 1:n
        rand() < p0 ? node_state[i] = I : node_state[i] = S
    end
    HyperNetwork(n, node_state)
end


# ====================================================================================
# ----------------------------- GRAPH MANIPULATION -----------------------------------


function add_hyperedge!(network::HyperNetwork, nodes::Vector{Int64})
    @assert all(nodes .<= get_num_nodes(network))
    
    vertices = Dict([(n, true) for n in nodes])
    SimpleHypergraphs.add_hyperedge!(network.hg; vertices = vertices)
    new_size = length(nodes)
    _add_to_hyperedge_dist!(network.hyperedge_dist, new_size)
    network.max_hyperedge_uid += 1
    push!(network.hyperedge_uid, network.max_hyperedge_uid)
end


function add_node!(network::HyperNetwork, hyperedges::Vector{Int64}, state::State)
    @assert all([h in network.hyperedge_uid for h in hyperedges])

    # update hyperedge_dist
    for h in hyperedges
        old_size = get_hyperedge_size(network, h)
        network.hyperedge_dist[old_size] -= 1
        new_size = old_size + 1
        _add_to_hyperedge_dist!(network.hyperedge_dist, new_size)
    end

    SH_hyperedges = Dict([(indexin(h, network.hyperedge_uid)[], true) for h in hyperedges])
    SimpleHypergraphs.add_vertex!(network.hg, hyperedges = SH_hyperedges, v_meta = state)
    network.state_dist[state] += 1
end


function _add_to_hyperedge_dist!(hyperedge_dist::Dict, new_size::Integer)
    if new_size in keys(hyperedge_dist)
        hyperedge_dist[new_size] += 1
    else
        hyperedge_dist[new_size] = 1
        # fill in all previous keys
        for size in new_size-1:-1:2
            if !(size in keys(hyperedge_dist))
                hyperedge_dist[size] = 0
            end
        end
    end
end


"""
remove_hyperedge!(network::HyperNetwork, hyperedge::Integer)
"""
function remove_hyperedge!(network::HyperNetwork, hyperedge::Integer)
    @assert hyperedge in network.hyperedge_uid
    num_hyperedges = get_num_hyperedges(network)
    
    # update hyperedge_dist
    old_size = get_hyperedge_size(network, hyperedge)
    network.hyperedge_dist[old_size] -= 1
    
    # update hyperedge_uid
    # The function SimpleHypergraphs.remove_hyperedge!() does not preserve the order 
    # of the hyperedges: when a hyperedge is deleted, the last column is moved to 
    # the posititon where the previous hyperedge was. 
    mid = indexin(hyperedge, network.hyperedge_uid)[]
    new_uid = pop!(network.hyperedge_uid)
    if mid != num_hyperedges
        network.hyperedge_uid[mid] = new_uid
    end
    
    SimpleHypergraphs.remove_hyperedge!(network.hg, mid)
    
    return nothing
end


"""
If the hyperedge is of size two, it is removed completely from the graph. 
"""
function remove_node_from_hyperedge!(network::HyperNetwork, node::Integer, hyperedge::Integer)
    old_size = get_hyperedge_size(network, hyperedge)
    if old_size == 2
        remove_hyperedge!(network, hyperedge)
    else
        network.hyperedge_dist[old_size] -= 1
        new_size = old_size - 1
        _add_to_hyperedge_dist!(network.hyperedge_dist, new_size)
        mid = indexin(hyperedge, network.hyperedge_uid)[]
        network.hg[node, mid] = nothing
    end
    return nothing
end


function set_state!(network::HyperNetwork, node::Integer, state::State)
    old_state = get_state(network, node)
    set_vertex_meta!(network.hg, state, node)
    network.state_dist[old_state] -= 1
    network.state_dist[state] += 1
end


function set_hyperedge_meta!(network::HyperNetwork, hyperedge::Integer, meta)
    mid = indexin(hyperedge, network.hyperedge_uid)[]
    SimpleHypergraphs.set_hyperedge_meta!(network.hg, meta, mid)
end


function set_hyperedge_meta!(network::HyperNetwork, hyperedge::Integer, key::Symbol, value::Any)
    meta = get_hyperedge_meta(network, hyperedge)
    meta[key] = value
    set_hyperedge_meta!(network, meta, hyperedge)
end


function get_hyperedge_meta(network::HyperNetwork, hyperedge::Integer)
    mid = indexin(hyperedge, network.hyperedge_uid)[]
    SimpleHypergraphs.get_hyperedge_meta(network.hg, mid)
end

function get_hyperedge_meta(network::HyperNetwork, hyperedge::Integer, key::Symbol)
    get_hyperedge_meta(network, hyperedge)[key]
end


# ====================================================================================
# --------------------------------- GRAPH INFO ---------------------------------------

function get_nodes(network::HyperNetwork)
    return 1:get_num_nodes(network)
end

function get_nodes(network::HyperNetwork, hyperedge::Integer)
    @assert hyperedge in network.hyperedge_uid
    mid = indexin(hyperedge, network.hyperedge_uid)[]
    return collect(keys(filter(d->d.second, getvertices(network.hg, mid))))
end

function get_hyperedges(network::HyperNetwork)
    return network.hyperedge_uid
end

function get_state(network::HyperNetwork, node::Integer)
    @assert node <= get_num_nodes(network)
    return get_vertex_meta(network.hg, node)
end

function get_node_to_state_dict(network::HyperNetwork)
    return Dict(node => get_state(network, node) for node in get_nodes(network))
end

function get_node_to_state_dict(network::HyperNetwork, hyperedge::Integer)
    @assert hyperedge in network.hyperedge_uid
    return Dict(node => get_state(network, node) for node in get_nodes(network, hyperedge))
end

function get_state_dist(network::HyperNetwork)
    return network.state_dist
end

function get_hyperedge_dist(network::HyperNetwork)
    return network.hyperedge_dist
end

function get_num_hyperedges(network::HyperNetwork)
    return nhe(network.hg)
end

function get_num_nodes(network::HyperNetwork)
    return nhv(network.hg)
end

function get_node_degree(network::HyperNetwork, node::Integer)
    @assert node <= get_num_nodes(network)
    return sum(values(network.hg.v2he[node]))
end

function get_hyperedge_size(network::HyperNetwork, hyperedge::Integer)
    @assert hyperedge in network.hyperedge_uid
    mid = indexin(hyperedge, network.hyperedge_uid)[]
    return sum(values(network.hg.he2v[mid]))
end

function get_max_hyperedge_size(network::HyperNetwork)
    return maximum([sum(values(d)) for d in network.hg.he2v])
end


"""
    is_active(network::HyperNetwork, hyperedge::Integer)

Returns true if the hyperedge contains nodes in different states, false if all states are equal. 
"""
function is_active(network::HyperNetwork, hyperedge::Integer)
    nodes = get_nodes(network, hyperedge)
    states = [get_state(network, n) for n in nodes]
    return length(unique(states)) > 1
end


"""
    get_twosection_graph(network::HyperNetwork)

Returns the two-section graph of the hypergraph as a SimpleGraph. 

A two-section of a hypergraph is a graph with the same vertices where two vertices 
are connected if they belong to the same hyperedge. Information about overlapping or 
parallel hyperedges is lost during conversion. 
"""
function get_twosection_graph(network::HyperNetwork)
    adjmatrix = get_twosection_adjacency_mx(network.hg, replace_weights=1)
    return Graphs.SimpleGraphs.SimpleGraph(adjmatrix)
end


# ====================================================================================
# ----------------------------- GRAPH CONSTRUCTION -----------------------------------


"""
    build_RSC_hg!(network::HyperNetwork, num_hyperedges::Tuple{Vararg{Integer}})

Populates the hypergraph with randomly distributed hyperedges of different dimensions. 

The number of hyperedges in each dimension is given by `num_hyperedges`, starting with d = 1 (e.g., edge between two nodes). 

The resulting hyperedges will contain only distinct nodes and each hyperedge occures only once.
It is however possible that some hyperedges will be subsets of others.
It is assumed that the hypergraph is empty; otherwise, the hyperedges will be added, but the conditions above are not guaranteed. 

The algorithm roughly follows the Iacopini paper, but uses the absolute number of hyperdeges instead of p_d and <k_d>.

TODO: make this more formal
"""
function build_RSC_hg!(network::HyperNetwork, num_hyperedges::Tuple{Vararg{Integer}})
    max_dim = length(num_hyperedges)
    n = get_num_nodes(network)
    for d in 1:max_dim
        num_inserted_hyperedges = 0
        # TODO: think of a better data structure if this ever becomes a bottleneck
        history::Vector{Vector{Int64}} = []
        while num_inserted_hyperedges < num_hyperedges[d]
            nodes = rand(1:n, d + 1)
            sort!(nodes)
            if allunique(nodes) && !(nodes in history)
                add_hyperedge!(network, nodes)
                push!(history, nodes)
                num_inserted_hyperedges += 1
            end
        end
    end
    return nothing
end


function build_RSC_hg_new!(network::HyperNetwork, num_hyperedges::Tuple{Vararg{Integer}})
    max_dim = length(num_hyperedges)
    n = get_num_nodes(network)
    for d in 1:max_dim
        # draw num_hyperedges[d] distinct indices of the combinations
        indices = rand(0:binomial(n - 1, d) - 1, num_hyperedges[d])
        for index in indices
            nodes = _index_to_combination(index, d)
            # _index_to_combination returns combinations with numbers starting from zero,
            # but we need them to start from one
            nodes .+= 1
            add_hyperedge!(network, nodes)
        end
    end
    return nothing
end


"""
Finds the combination of size `size` at the given index. 

Each combination is a set of unique numbers greater or equal than zero sorted in *decreasing* order. 
If all combinations are sorted in a lexographic order, a unique index can be assigned to every combination. 
This function computes the reverse mapping: given an index, it finds the corresponding combination. 

See https://en.wikipedia.org/wiki/Combinatorial_number_system#Finding_the_k-combination_for_a_given_number
"""
function _index_to_combination(index::Integer, size::Integer)
    @assert index ≥ binomial(size - 1, size)
    combination = []
    for k = size:-1:1
        ck = k - 1
        binom = binomial(ck + 1, k)
        while binom <= index
            ck += 1
            binom *= (ck + 1) 
            binom ÷= (ck + 1 - k)
        end
        push!(combination, ck)
        index -= binomial(ck, k)
    end
    return combination
end


"""
    build_regular_hg!(network::HyperNetwork, degrees::Tuple{Vararg{Integer}})

Populates the hypergraph with hyperedges such that every node has degrees {d_1, d_2, ...}. 
"""
function build_regular_hg!(network::HyperNetwork, degrees::Tuple{Vararg{Integer}})
# TODO
end
