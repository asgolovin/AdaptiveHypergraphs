using SimpleHypergraphs
using StatsBase
using Random

@enum State::Bool S=false I=true

"""
    HyperNetwork

A higher-order network where every node can be in a particular state (SIS or SIR).

The topology of the network is represented by a hypergraph: a graph, where every hyperedge can connect multiple vertices, and not just two. 
"""
mutable struct HyperNetwork
    # The underlying hypergraph
    hg::Hypergraph{Bool, State}
    # Number of nodes in a particular state
    state_dist::Dict{State, Integer}
end

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
    println(node_state)
    HyperNetwork(n, node_state)
end


function add_hyperedge!(network::HyperNetwork, vertices)
    SimpleHypergraphs.add_hyperedge!(network.hg; vertices = Dict([(v, true) for v in vertices]))
end

function add_node!(network::HyperNetwork, hyperedges, state::State)
    SimpleHypergraphs.add_vertex!(network.hg, hyperedges = Dict([(h, true) for h in hyperedges]), v_meta = state)
    network.state_dist[state] += 1
end

function remove_hyperedge!(network::HyperNetwork, hyperedge)
    @assert hyperedge <= network.hg.nhe()
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