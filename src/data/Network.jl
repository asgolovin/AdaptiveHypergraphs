using SimpleHypergraphs
using StatsBase
using Random
using Graphs
using Combinatorics

export HyperNetwork,
       add_hyperedge!, include_node!, delete_hyperedge!, remove_node!, set_state!,
       get_nodes, get_hyperedges, get_intersecting_hyperedges, get_state, get_state_map,
       get_state_count, get_hyperedge_dist, get_num_hyperedges, get_num_active_hyperedges,
       get_num_nodes, get_node_degree, get_hyperedge_size, get_max_size, get_motif_count,
       is_active, get_twosection_graph, build_regular_hg!, build_RSC_hg!

"""
A node is included multiple times into a hyperdege.
"""
struct DegenerateHyperedge <: Exception
    node::Int64
    hyperedge::Int64
end

function Base.showerror(io::IO, e::DegenerateHyperedge)
    return print(io,
                 "Node $(e.node) already belongs to the hyperdege $(e.hyperedge) and cannot be added again.")
end

"""
A hyperedge with this set of nodes exists already.
"""
struct ParallelHyperedge <: Exception end

"""
    HyperNetwork

A higher-order network where every node can be in a particular state.

The topology of the network is represented by a hypergraph: a graph, where every hyperedge can connect multiple vertices, and not just two. 

Other than the hypergraph and the states, the struct keeps track of any statistics related to the system such as the number of infected nodes. 

Important assumptions / modifications to traditional hypergraphs:
 - Empty hyperedges or hyperedges with just one node can't exist. If a node is removed from a hyperedge of size two, the hyperedge is deleted.
 - The number of nodes stays constant. Even if a node is not connected to any hyperedeges, it is not deleted. 
 - Degenerate hyperedges (hyperedges where the same node occures multiple times) can't exist

Assumptions which are not yet quite enforced in the code, but are a TODO:
 - Parallel hyperedges (hyperedges with exactly the same set of nodes) can't exist
"""
mutable struct HyperNetwork
    # The underlying hypergraph
    hg::Hypergraph{Bool,State}
    # Number of nodes
    num_nodes::Int64
    # maximum allowed size of hyperedges 
    max_size::Int64
    # whether to keep track of the number of motifs (expensive!)
    track_motif_count::Bool
    # Motif -> number of motifs
    motif_count::Dict{AbstractMotif,Int64}
    # Size -> number of active hyperedges
    active_hyperedges::Dict{Int64,Int64}
    # Number of nodes in a particular state
    state_count::Dict{State,Int64}
    # Size -> number of hyperdeges
    hyperedge_dist::Dict{Int64,Int64}
    # A map of UIDs to the size of the hyperedege
    hyperedge_size::Dict{Int64,Int64}
    # A bijective map where the vector indices correspond to the indices of the columns 
    # in the incidence matrix of the hypergraph (matrix IDs or MIDs) and the values 
    # to unique IDs (UIDs). The UIDs start with 1. 
    # Only the UIDs are exposed outside of the class; the MIDs should 
    # only be used to communicate with SimpleHypergraphs. 
    # This mapping is important because the native hyperedge indices from 
    # SimpleHypergraphs are not preserved when a hyperedge is deleted. 
    hyperedge_uid::Vector{Int64}
    # A reverse map of UIDs to MIDs.
    uid_to_mid::Dict{Int64,Int64}
    # current highest hyperedge UID
    max_hyperedge_uid::Int64
end

# ====================================================================================
# ------------------------------- CONSTRUCTORS ---------------------------------------

"""
    HyperNetwork(n::Int64, node_state::Vector{Union{Nothing, State}}, max_size::Int64)

Create an empty network with `n` nodes and no hyperedges.
`node_state` denotes the state of each node. 
`max_size` is the maximum allowed hyperedge size. 

The optional flag `track_motif_count` allows to switch on expensive tracking of the number of motifs.  
"""
function HyperNetwork(n::Int64,
                      node_state::Vector{Union{Nothing,State}},
                      max_size::Int64;
                      track_motif_count::Bool=false)
    @assert length(node_state) == n
    @assert max_size >= 2

    # create the hypergraph
    matrix = Matrix{Union{Nothing,Bool}}(nothing, (n, 0))
    hg = Hypergraph{Bool,State}(matrix; v_meta=node_state)

    # initialize empty data structures
    motif_count = Dict{AbstractMotif,Int64}([l => 0 for l in all_motifs(max_size)])
    active_hyperedges = Dict([size => 0 for size in 2:max_size])
    state_count = countmap(node_state)
    # fill the missing keys
    for state in instances(State)
        if state ∉ keys(state_count)
            state_count[state] = 0
        end
    end
    motif_count[OrderZeroMotif(A)] = state_count[A]
    motif_count[OrderZeroMotif(B)] = state_count[B]
    hyperedge_dist = Dict([size => 0 for size in 2:max_size])
    hyperedge_size = Dict{Int64,Int64}()
    hyperedge_uid = Vector{Int64}()
    uid_to_mid = Dict{Int64,Int64}()

    return HyperNetwork(hg, n, max_size, track_motif_count,
                        motif_count, active_hyperedges, state_count,
                        hyperedge_dist, hyperedge_size,
                        hyperedge_uid, uid_to_mid, 0)
end

"""
    HyperNetwork(n::Int64)

Create an empty network with `n` nodes and no hyperedges.
All nodes have the opinion A.
"""
function HyperNetwork(n::Int64, max_size::Int64; track_motif_count::Bool=false)
    node_state = Vector{Union{Nothing,State}}(nothing, n)
    fill!(node_state, A)
    return HyperNetwork(n, node_state, max_size; track_motif_count)
end

"""
    HyperNetwork(n::Int64, p0::Float64)

Create an empty network with n nodes and no hyperedges.
Each node has opinion A with probability p0. 
"""
function HyperNetwork(n::Int64, p0::Float64, max_size::Int64; track_motif_count::Bool=false)
    node_state = Vector{Union{Nothing,State}}(nothing, n)
    for i in 1:n
        rand() < p0 ? node_state[i] = A : node_state[i] = B
    end
    return HyperNetwork(n, node_state, max_size; track_motif_count)
end

function Base.show(io::IO, network::HyperNetwork)
    num_nodes = get_num_nodes(network)
    return print(io, "HyperNetwork($num_nodes)")
end

function Base.show(io::IO, ::MIME"text/plain", network::HyperNetwork)
    num_nodes = get_num_nodes(network)
    num_hyperedges = get_num_hyperedges(network)

    println(io, "HyperNetwork with $num_nodes nodes and $num_hyperedges hyperedges")
    println(io, "states:")
    for state in instances(State)
        println(io, "  $state => $(network.state_count[state]) nodes")
    end
    println(io, "hyperedges:")
    for size in 2:get_max_size(network)
        num_hyperedges = network.hyperedge_dist[size]
        num_active = network.active_hyperedges[size]
        println(io,
                "  size $size => $num_hyperedges hyperdeges, $num_active/$num_hyperedges active")
    end
    println(io, "motifs:")
    for motif in all_motifs(get_max_size(network))
        num_motifs = network.motif_count[motif]
        println(io, "  $motif => $num_motifs motifs")
    end
end

# ====================================================================================
# ----------------------------- GRAPH MANIPULATION -----------------------------------

"""
    add_hyperedge!(network::HyperNetwork, nodes)

Create a new hyperedge with nodes `nodes` and add it to `network`.
"""
function add_hyperedge!(network::HyperNetwork, nodes)
    @assert all(nodes .<= get_num_nodes(network))
    @assert length(nodes) >= 2
    @assert allunique(nodes)
    @assert length(nodes) <= network.max_size

    # add the hyperedge to the hypergraph
    vertices = Dict([(n, true) for n in nodes])
    mid = SimpleHypergraphs.add_hyperedge!(network.hg; vertices=vertices)

    # update hyperedge_dist
    new_size = length(nodes)
    network.hyperedge_dist[new_size] += 1

    # update uid
    network.max_hyperedge_uid += 1
    uid = network.max_hyperedge_uid
    push!(network.hyperedge_uid, uid)
    network.uid_to_mid[uid] = mid
    network.hyperedge_size[uid] = length(nodes)
    hyperedge = uid

    # update active_hyperedeges
    if is_active(network, hyperedge)
        network.active_hyperedges[new_size] += 1
    end

    # increase new values in motif count
    if network.track_motif_count
        motif = _get_1st_order_motif(network, hyperedge)
        network.motif_count[motif] += 1
        for motif in _get_2nd_order_motifs(network, hyperedge)
            network.motif_count[motif] += 1
        end
    end

    return network.max_hyperedge_uid
end

"""
    include_node!(network::HyperNetwork, node::Int64, hyperedge::Int64)

Add an existing node to an existing hyperedge. 
"""
function include_node!(network::HyperNetwork, node::Int64, hyperedge::Int64)
    @assert hyperedge in network.hyperedge_uid
    @assert 1 <= node <= get_num_nodes(network)
    @assert length(get_nodes(network, hyperedge)) + 1 <= network.max_size
    if node in get_nodes(network, hyperedge)
        throw(DegenerateHyperedge(node, hyperedge))
    end

    # decrease old values in motif_count
    if network.track_motif_count
        motif = _get_1st_order_motif(network, hyperedge)
        network.motif_count[motif] -= 1
        for motif in _get_2nd_order_motifs(network, hyperedge)
            network.motif_count[motif] -= 1
        end
    end

    # update hyperedge_dist
    old_size = get_hyperedge_size(network, hyperedge)
    network.hyperedge_dist[old_size] -= 1
    new_size = old_size + 1
    network.hyperedge_dist[new_size] += 1
    network.hyperedge_size[hyperedge] += 1

    active_before = is_active(network, hyperedge)

    # update the network
    mid = network.uid_to_mid[hyperedge]
    network.hg[node, mid] = true

    # increase new values in motif_count
    if network.track_motif_count
        motif = _get_1st_order_motif(network, hyperedge)
        network.motif_count[motif] += 1
        for motif in _get_2nd_order_motifs(network, hyperedge)
            network.motif_count[motif] += 1
        end
    end

    # update active_hyperedeges if the hyperedge is now active
    if is_active(network, hyperedge)
        if active_before # if it was already active
            network.active_hyperedges[new_size] += 1
            network.active_hyperedges[old_size] -= 1
        else # if it became active
            network.active_hyperedges[new_size] += 1
        end
    end

    return network
end

"""
    delete_hyperedge!(network::HyperNetwork, hyperedge::Int64)

Remove `hyperedge` from ̀`network`.
"""
function delete_hyperedge!(network::HyperNetwork, hyperedge::Int64)
    @assert hyperedge in network.hyperedge_uid
    num_hyperedges = get_num_hyperedges(network)

    # decrease old values in motif_count
    if network.track_motif_count
        motif = _get_1st_order_motif(network, hyperedge)
        network.motif_count[motif] -= 1
        for motif in _get_2nd_order_motifs(network, hyperedge)
            network.motif_count[motif] -= 1
        end
    end

    # update hyperedge_dist
    old_size = get_hyperedge_size(network, hyperedge)
    network.hyperedge_dist[old_size] -= 1

    # update active_hyperedeges if the hyperedge was active
    if is_active(network, hyperedge)
        network.active_hyperedges[old_size] -= 1
    end

    # update hyperedge_uid
    # The function SimpleHypergraphs.remove_hyperedge!() does not preserve the order 
    # of the hyperedges: when a hyperedge is deleted, the last column is moved to 
    # the posititon where the previous hyperedge was. 
    mid = network.uid_to_mid[hyperedge]
    delete!(network.uid_to_mid, hyperedge)
    delete!(network.hyperedge_size, hyperedge)
    new_uid = pop!(network.hyperedge_uid)
    if mid != num_hyperedges
        network.hyperedge_uid[mid] = new_uid
        network.uid_to_mid[new_uid] = mid
    end

    SimpleHypergraphs.remove_hyperedge!(network.hg, mid)

    return network
end

"""
    remove_node!(network::HyperNetwork, node::Int64,
                        hyperedge::Int64)

Remove `node` from `hyperedge`.

If the hyperedge was of size two, it is deleted completely from the graph (hyperdeges of size one 
are not allowed). However, the node continues to exist even if it is not attached to any 
hyperedeges anymore. 
"""
function remove_node!(network::HyperNetwork, node::Int64,
                      hyperedge::Int64)
    old_size = get_hyperedge_size(network, hyperedge)
    if old_size == 2
        delete_hyperedge!(network, hyperedge)
        return network
    end

    # decrease old values in motif_count
    if network.track_motif_count
        motif = _get_1st_order_motif(network, hyperedge)
        network.motif_count[motif] -= 1
        for motif in _get_2nd_order_motifs(network, hyperedge)
            network.motif_count[motif] -= 1
        end
    end

    active_before = is_active(network, hyperedge)
    network.hyperedge_dist[old_size] -= 1
    network.hyperedge_size[hyperedge] -= 1
    new_size = old_size - 1
    network.hyperedge_dist[new_size] += 1
    mid = network.uid_to_mid[hyperedge]
    network.hg[node, mid] = nothing

    # increase new values in motif_count
    if network.track_motif_count
        motif = _get_1st_order_motif(network, hyperedge)
        network.motif_count[motif] += 1
        for motif in _get_2nd_order_motifs(network, hyperedge)
            network.motif_count[motif] += 1
        end
    end

    # update active_hyperedeges if it was active before
    if active_before
        if is_active(network, hyperedge) # stayed active
            network.active_hyperedges[old_size] -= 1
            network.active_hyperedges[new_size] += 1
        else # switched from active to inactive
            network.active_hyperedges[old_size] -= 1
        end
    end
    return network
end

"""
    set_state!(network::HyperNetwork, node::Int64, state::State)

Set the state of `node` to `state`.
"""
function set_state!(network::HyperNetwork, node::Int64, state::State)
    @assert 1 <= node <= get_num_nodes(network)

    hyperedges = get_hyperedges(network, node)
    active_before = [is_active(network, h) for h in hyperedges]
    old_state = get_state(network, node)

    # decrease old values in motif_count
    if network.track_motif_count
        for hyperedge in get_hyperedges(network, node)
            motif = _get_1st_order_motif(network, hyperedge)
            network.motif_count[motif] -= 1
        end
        for motif in _get_2nd_order_node_motifs(network, node)
            network.motif_count[motif] -= 1
        end
    end

    set_vertex_meta!(network.hg, state, node)

    network.state_count[old_state] -= 1
    network.state_count[state] += 1
    network.motif_count[OrderZeroMotif(A)] = network.state_count[A]
    network.motif_count[OrderZeroMotif(B)] = network.state_count[B]

    # increase old values in motif_count
    if network.track_motif_count
        for hyperedge in get_hyperedges(network, node)
            motif = _get_1st_order_motif(network, hyperedge)
            network.motif_count[motif] += 1
        end
        for motif in _get_2nd_order_node_motifs(network, node)
            network.motif_count[motif] += 1
        end
    end

    for (i, h) in enumerate(hyperedges)
        size = get_hyperedge_size(network, h)
        if !active_before[i] && is_active(network, h)
            # if the edge switched from _inactive to active_
            network.active_hyperedges[size] += 1

        elseif active_before[i] && !is_active(network, h)
            # if the edge switched from _active to inactive_
            network.active_hyperedges[size] -= 1
        end
    end

    return network
end

# ====================================================================================
# --------------------------------- GRAPH INFO ---------------------------------------

function get_nodes(network::HyperNetwork)
    return collect(1:(network.num_nodes))
end

function get_nodes(network::HyperNetwork, hyperedge::Int64)
    @assert hyperedge in network.hyperedge_uid
    mid = network.uid_to_mid[hyperedge]
    return collect(keys(filter(d -> d.second, getvertices(network.hg, mid))))
end

function get_hyperedges(network::HyperNetwork)
    return copy(network.hyperedge_uid)
end

function get_hyperedges(network::HyperNetwork, node::Int64)
    @assert 1 <= node <= get_num_nodes(network)
    mids = collect(keys(filter(d -> d.second, gethyperedges(network.hg, node))))
    return network.hyperedge_uid[mids]
end

"""
Return all hyperedges which intersect with the given `hyperedge`
"""
function get_intersecting_hyperedges(network::HyperNetwork, hyperedge::Int64)
    @assert hyperedge in network.hyperedge_uid
    intersecting_hyperedges = []
    for node in get_nodes(network, hyperedge)
        append!(intersecting_hyperedges, get_hyperedges(network, node))
    end
    unique!(intersecting_hyperedges)
    # remove the hyperedge itself from the list
    filter!(x -> x != hyperedge, intersecting_hyperedges)
    return intersecting_hyperedges
end

"""
Return the first-order motif corresponding to the hyperedge `hyperedge`
"""
function _get_1st_order_motif(network::HyperNetwork, hyperedge::Int64)
    statecount = get_state_count(network, hyperedge)
    return OrderOneMotif(statecount[A], statecount[B])
end

"""
Return the second-order motif formed by two intersecting hyperedges `hyperedge1` and `hyperedge2`.
"""
function _get_2nd_order_motif(network::HyperNetwork, hyperedge1::Int64, hyperedge2::Int64)
    statecount1 = get_state_count(network, hyperedge1)
    statecount2 = get_state_count(network, hyperedge2)
    # collect all nodes in which the hyperedges intersect
    int_nodes = []
    for node in get_nodes(network, hyperedge1)
        if node ∈ get_nodes(network, hyperedge2)
            push!(int_nodes, node)
        end
    end
    int_statecount = countmap([get_state(network, node) for node in int_nodes])
    for state in instances(State)
        if state ∉ keys(int_statecount)
            int_statecount[state] = 0
        end
    end
    statecount1[A] -= int_statecount[A]
    statecount2[A] -= int_statecount[A]
    statecount1[B] -= int_statecount[B]
    statecount2[B] -= int_statecount[B]

    motif = OrderTwoMotif((statecount1[A], statecount1[B]),
                          (int_statecount[A], int_statecount[B]),
                          (statecount2[A], statecount2[B]))

    return motif
end

"""
Return all second-order motifs which contain the hyperedge `hyperedge`. 
"""
function _get_2nd_order_motifs(network::HyperNetwork, hyperedge::Int64)
    motifs = []
    for intersecting_he in get_intersecting_hyperedges(network, hyperedge)
        motif = _get_2nd_order_motif(network, hyperedge, intersecting_he)
        push!(motifs, motif)
    end
    return motifs
end

"""
Return all second-order motifs which contain the given node.
"""
function _get_2nd_order_node_motifs(network::HyperNetwork, node::Int64)
    motifs = []
    # hyperedges which include the node
    hyperedges = get_hyperedges(network, node)

    # add all triples which intersect in the node
    int_state = get_state(network, node)
    for (i, he1) in enumerate(hyperedges)
        for he2 in hyperedges[(i + 1):end]
            motif = _get_2nd_order_motif(network, he1, he2)
            push!(motifs, motif)
        end
    end

    # add all triples which *do not* intersect in the node
    for he1 in hyperedges
        for he2 in get_intersecting_hyperedges(network, he1)
            if node ∈ get_nodes(network, he2)
                continue
            end
            motif = _get_2nd_order_motif(network, he1, he2)
            push!(motifs, motif)
        end
    end
    return motifs
end

function get_state(network::HyperNetwork, node::Int64)
    @assert 1 <= node <= get_num_nodes(network)
    return SimpleHypergraphs.get_vertex_meta(network.hg, node)
end

"""
    get_state_map(network::HyperNetwork)

Return a dict which maps every node to its state.
"""
function get_state_map(network::HyperNetwork)
    return Dict(node => get_state(network, node) for node in get_nodes(network))
end

"""
    get_state_map(network::HyperNetwork, hyperedge::Int64)
    
Return a dict which maps every node in `hyperedge` to its state.
"""
function get_state_map(network::HyperNetwork, hyperedge::Int64)
    @assert hyperedge in network.hyperedge_uid
    statemap = Dict(node => get_state(network, node)
                    for node in get_nodes(network, hyperedge))
    return statemap
end

"""
    get_state_count(network::HyperNetwork)

Return a dict which maps every state to the number of nodes in this state. 
"""
function get_state_count(network::HyperNetwork)
    return copy(network.state_count)
end

"""
    get_state_count(network::HyperNetwork)

Return a dict which maps every state of every node in a hyperdege to the number of nodes in this state. 
"""
function get_state_count(network::HyperNetwork, hyperedge::Int64)
    nodes = get_nodes(network, hyperedge)
    map = countmap([get_state(network, node) for node in nodes])
    for state in instances(State)
        if state ∉ keys(map)
            map[state] = 0
        end
    end
    return map
end

"""
    get_hyperedge_dist(network::HyperNetwork)

Return a dict which maps the sizes of hyperdeges to the number of hyperedeges of this size in the graph.
"""
function get_hyperedge_dist(network::HyperNetwork)
    return copy(network.hyperedge_dist)
end

function get_num_hyperedges(network::HyperNetwork)
    return nhe(network.hg)
end

function get_num_active_hyperedges(network::HyperNetwork)
    return sum(values(network.active_hyperedges))
end

function get_num_active_hyperedges(network::HyperNetwork, size::Int64)
    return network.active_hyperedges[size]
end

function get_num_nodes(network::HyperNetwork)
    return network.num_nodes
end

function get_node_degree(network::HyperNetwork, node::Int64)
    @assert 1 <= node <= get_num_nodes(network)
    return sum(values(network.hg.v2he[node]))
end

function get_hyperedge_size(network::HyperNetwork, hyperedge::Int64)
    # @assert hyperedge in network.hyperedge_uid
    return network.hyperedge_size[hyperedge]
end

"""
    get_max_size(network::HyperNetwork)

Return the maximum *historical* hyperedge size. 

The "historical" part is important for example for functions which plot the evolution
of hyperedge sizes over time. 
"""
function get_max_size(network::HyperNetwork)
    return maximum(keys(network.hyperedge_dist))
end

"""
    get_motif_count(network::HyperNetwork)

Return a dict mapping motifs to the number of occurances.
"""
function get_motif_count(network::HyperNetwork)
    if !network.track_motif_count
        throw(ErrorException("The motif count is not tracked. Please initialize the network with track_motif_count = true."))
    end
    return network.motif_count
end

"""
    is_active(network::HyperNetwork, hyperedge::Int64)

Return true if the hyperedge contains nodes in different states, false if all states are equal. 
"""
function is_active(network::HyperNetwork, hyperedge::Int64)
    nodes = get_nodes(network, hyperedge)
    states = [get_state(network, n) for n in nodes]
    return length(unique(states)) > 1
end

"""
    get_twosection_graph(network::HyperNetwork)

Return the two-section graph of the hypergraph as a SimpleGraph. 

A two-section of a hypergraph is a graph with the same vertices where two vertices 
are connected if they belong to the same hyperedge. Information about overlapping or 
parallel hyperedges is lost during conversion. 
"""
function get_twosection_graph(network::HyperNetwork)
    adjmatrix = get_twosection_adjacency_mx(network.hg; replace_weights=1)
    return Graphs.SimpleGraphs.SimpleGraph(adjmatrix)
end

# ====================================================================================
# ----------------------------- GRAPH CONSTRUCTION -----------------------------------

"""
    build_RSC_hg!(network::HyperNetwork, num_hyperedges::Tuple{Vararg{Int64}})

Populate the hypergraph with randomly distributed hyperedges of different dimensions. 

The number of hyperedges in each dimension is given by `num_hyperedges`, starting with d = 1 (e.g., edge between two nodes). 

The resulting hyperedges will contain only distinct nodes and each hyperedge occures only once.
It is however possible that some hyperedges will be subsets of others.
It is assumed that the hypergraph is empty; otherwise, the hyperedges will be added, but the conditions above are not guaranteed. 

The algorithm roughly follows the Iacopini paper, but uses the absolute number of hyperdeges instead of p_d and <k_d>.
"""
function build_RSC_hg!(network::HyperNetwork, num_hyperedges::Tuple{Vararg{Int64}})
    max_dim = length(num_hyperedges)
    n = get_num_nodes(network)
    for size in 2:get_max_size(network)
        network.active_hyperedges[size] = 0
        network.hyperedge_dist[size] = 0
    end
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
    return network
end

"""
A fancy variant of `build_RSC_hg!` using the combinatorial number system. 
The effect is the same, only the algorithm is different.

Slower than `build_RSC_hg!`. 
"""
function build_RSC_hg_new!(network::HyperNetwork, num_hyperedges::Tuple{Vararg{Int64}})
    max_dim = length(num_hyperedges)
    n = get_num_nodes(network)
    for d in 1:max_dim
        # draw num_hyperedges[d] distinct indices of the combinations
        indices = rand(0:(binomial(n - 1, d) - 1), num_hyperedges[d])
        for index in indices
            nodes = _index_to_combination(index, d)
            # _index_to_combination returns combinations with numbers starting from zero,
            # but we need them to start from one
            nodes .+= 1
            add_hyperedge!(network, nodes)
        end
    end
    return network
end

"""
Find the combination of size `size` at the given index. 

Each combination is a set of unique numbers greater or equal than zero sorted in *decreasing* order. 
If all combinations are sorted in a lexographic order, a unique index can be assigned to every combination. 
This function computes the reverse mapping: given an index, it finds the corresponding combination. 

See https://en.wikipedia.org/wiki/Combinatorial_number_system#Finding_the_k-combination_for_a_given_number
"""
function _index_to_combination(index::Int64, size::Int64)
    @assert index ≥ binomial(size - 1, size)
    combination = []
    for k in size:-1:1
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