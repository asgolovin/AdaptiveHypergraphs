using SimpleHypergraphs
using StatsBase
using Random

@enum State::Bool S=false I=true

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
end


# ====================================================================================
# ------------------------------- CONSTRUCTORS ---------------------------------------


"""
    HyperNetwork(n::Integer, node_state::Vector{Union{Nothing, State}})

Create an empty network with n nodes and no hyperedges.
`node_state` denotes the state of each node. 
"""
function HyperNetwork(n::Integer,
                      node_state::Vector{Union{Nothing, State}}
                      )
    @assert length(node_state) == n
    matrix = Matrix{Union{Nothing, Bool}}(nothing, (n, 0))
    hg = Hypergraph{Bool, State}(matrix; v_meta=node_state)
    state_dist = countmap(node_state)
    HyperNetwork(hg, state_dist)
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
    HyperNetwork(hg, Dict(S => n, I => 0))
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


function add_hyperedge!(network::HyperNetwork, nodes)
    @assert all(nodes .<= get_num_nodes(network))
    SimpleHypergraphs.add_hyperedge!(network.hg; vertices = Dict([(n, true) for n in nodes]))
end


function add_node!(network::HyperNetwork, hyperedges, state::State)
    @assert all(hyperedges .<= get_num_hyperedges(network)) 
    SimpleHypergraphs.add_vertex!(network.hg, hyperedges = Dict([(h, true) for h in hyperedges]), v_meta = state)
    network.state_dist[state] += 1
end


function remove_hyperedge!(network::HyperNetwork, hyperedge)
    @assert hyperedge <= get_num_nodes(network)
    SimpleHypergraphs.remove_hyperedge!(network.hg, hyperedge)
end


function set_state!(network::HyperNetwork, node, state::State)
    old_state = get_state(network, node)
    set_vertex_meta!(network.hg, state, node)
    network.state_dist[old_state] -= 1
    network.state_dist[state] += 1
end


function get_state(network::HyperNetwork, node)
    return get_vertex_meta(network.hg, node)
end


function get_num_hyperedges(network::HyperNetwork)
    return nhe(network.hg)
end


function get_num_nodes(network::HyperNetwork)
    return nhv(network.hg)
end


function get_node_degree(network::HyperNetwork, node::Integer)
    return length(gethyperedges(network.hg, node))
end

function get_nodes(network::HyperNetwork, hyperedge::Integer)
    return get_vertices(network.hg, hyperedge)
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
            # TODO: maybe move this to add_hyperedge? 
            if allunique(nodes) && !(nodes in history)
                add_hyperedge!(network, nodes)
                push!(history, nodes)
                num_inserted_hyperedges += 1
            end
        end
    end
    return nothing
end


"""
    build_regular_hg!(network::HyperNetwork, degrees::Tuple{Vararg{Integer}})

Populates the hypergraph with hyperedges such that every node has degrees {d_1, d_2, ...}. 
"""
function build_regular_hg!(network::HyperNetwork, degrees::Tuple{Vararg{Integer}})
# TODO
end
