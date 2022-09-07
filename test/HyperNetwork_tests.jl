using Random
using StatsBase

@testset "HyperNetwork: constructors" begin
    # create an empty network with all suseptible nodes
    n = 5
    network = HyperNetwork(n)
    @test get_num_nodes(network) == n
    @test get_num_hyperedges(network) == 0

    # an empty network with a given distribution of infected nodes
    node_state = Vector{Union{Nothing,State}}(nothing, n)
    fill!(node_state, AdaptiveHypergraphs.S)
    node_state[2] = AdaptiveHypergraphs.I
    node_state[5] = AdaptiveHypergraphs.I
    network = HyperNetwork(n, node_state)
    @test get_state_count(network)[AdaptiveHypergraphs.S] == n - 2
    @test get_state_count(network)[AdaptiveHypergraphs.I] == 2

    # an empty network with a certain probability of infection
    Random.seed!(42)
    network = HyperNetwork(n, 0.5)
    @test get_state_count(network)[AdaptiveHypergraphs.S] == 3
    @test get_state_count(network)[AdaptiveHypergraphs.I] == 2
end

@testset "Hypernetwork: graph manipulation" begin
    n = 5
    node_state = Vector{Union{Nothing,State}}(nothing, n)
    fill!(node_state, AdaptiveHypergraphs.S)
    node_state[1] = AdaptiveHypergraphs.I
    node_state[3] = AdaptiveHypergraphs.I
    network = HyperNetwork(n, node_state)

    add_hyperedge!(network, (1, 2))
    add_hyperedge!(network, (2, 3))
    add_hyperedge!(network, (1, 3, 4))
    # hyperedge 1: [1, 2], hyperedge 2: [2, 3], hyperedge 3: [1, 3, 4]

    @test get_num_nodes(network) == n
    @test get_num_hyperedges(network) == 3
    @test get_node_degree(network, 1) == 2
    @test get_node_degree(network, 5) == 0
    hyperedge_dist = get_hyperedge_dist(network)
    @test hyperedge_dist[2] == 2
    @test hyperedge_dist[3] == 1

    @test_throws AssertionError add_hyperedge!(network, (1,))
    @test_throws AssertionError add_hyperedge!(network, (1, 1))
    @test_throws AssertionError add_hyperedge!(network, (42, 234))
    @test_throws AssertionError delete_hyperedge!(network, 213)

    include_node!(network, 3, 1)
    include_node!(network, 5, 3)
    # hyperedge 1: [1, 2, 3], hyperedge 2: [2, 3], hyperedge 3: [1, 3, 4, 5]

    @test all([1, 2, 3] .∈ Ref(get_nodes(network, 1)))
    @test all([1, 3, 4, 5] .∈ Ref(get_nodes(network, 3)))
    hyperedge_dist = get_hyperedge_dist(network)
    @test hyperedge_dist[2] == 1
    @test hyperedge_dist[3] == 1
    @test hyperedge_dist[4] == 1

    delete_hyperedge!(network, 3)
    # hyperedge 1: [1, 2, 3], hyperedge 2: [2, 3]
    @test get_num_hyperedges(network) == 2
    hyperedge_dist = get_hyperedge_dist(network)
    @test hyperedge_dist[4] == 0
    # get_max_hyperedge_size still returns the historical max value
    @test get_max_hyperedge_size(network) == 4

    # removing a node from a hyperedge of size 2 deletes this hyperedge
    remove_node!(network, 2, 2)
    # hyperedge 1: [1, 2, 3]
    @test get_num_hyperedges(network) == 1
    # can't get nodes of a non-existing hyperedge
    @test_throws AssertionError get_nodes(network, 2)
    @test get_node_degree(network, 3) == 1
    @test get_hyperedges(network, 3) == [1]

    # removing a node from a larger hyperedge preserves the hyperedge
    remove_node!(network, 1, 1)
    # hyperedge 1: [2, 3]
    @test get_hyperedge_size(network, 1) == 2
    @test get_node_degree(network, 1) == 0
    @test get_hyperedges(network, 1) == []

    # check that after all manipulations, the number of active hyperedges is still correct
    @test get_num_active_hyperedges(network) == 1
end

@testset "Hypernetwork: graph info" begin
    n = 5
    network = HyperNetwork(n)
    add_hyperedge!(network, (1, 3, 4))

    # test that all nodes are in the same state S (i.e., the hyperedge is not active)
    @test is_active(network, 1) == false
    @test get_num_active_hyperedges(network) == 0
    @test get_state_count(network) ==
          Dict(AdaptiveHypergraphs.S => 5, AdaptiveHypergraphs.I => 0)

    set_state!(network, 3, AdaptiveHypergraphs.I)
    # now this is no longer the case
    @test is_active(network, 1) == true
    @test get_num_active_hyperedges(network) == 1

    @test get_state_map(network) ==
          Dict(1 => AdaptiveHypergraphs.S, 2 => AdaptiveHypergraphs.S,
               3 => AdaptiveHypergraphs.I, 4 => AdaptiveHypergraphs.S,
               5 => AdaptiveHypergraphs.S)
    @test get_state_map(network, 1) ==
          Dict(1 => AdaptiveHypergraphs.S, 3 => AdaptiveHypergraphs.I,
               4 => AdaptiveHypergraphs.S)
    @test get_state_count(network) ==
          Dict(AdaptiveHypergraphs.S => 4, AdaptiveHypergraphs.I => 1)
end

@testset "Hypernetwork: graph construction" begin
    n = 50
    network = HyperNetwork(n, 0.4)
    build_RSC_hg!(network, (10, 20, 30))

    @test get_num_hyperedges(network) == 60

    # check that the number of active hyperedges is correct
    num_active = 0
    hyperedge_dist = Dict(2 => 0, 3 => 0, 4 => 0)
    hyperedge_size = Dict()
    for h in get_hyperedges(network)
        if is_active(network, h)
            num_active += 1
        end
        hyperedge_dist[get_hyperedge_size(network, h)] += 1
        hyperedge_size[h] = get_hyperedge_size(network, h)
    end
    @test num_active == get_num_active_hyperedges(network)
    @test hyperedge_dist == get_hyperedge_dist(network)
    @test countmap(values(hyperedge_size)) == hyperedge_dist
end