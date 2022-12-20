using Random
using StatsBase

@testset "HyperNetwork: constructors" begin
    # create an empty network with all suseptible nodes
    n = 5
    network = HyperNetwork(n, 2)
    @test get_num_nodes(network) == n
    @test get_num_hyperedges(network) == 0

    # an empty network with a given distribution of infected nodes
    node_state = Vector{Union{Nothing,State}}(nothing, n)
    fill!(node_state, AdaptiveHypergraphs.A)
    node_state[2] = AdaptiveHypergraphs.B
    node_state[5] = AdaptiveHypergraphs.B
    network = HyperNetwork(n, node_state, 2)
    @test get_state_count(network)[AdaptiveHypergraphs.A] == n - 2
    @test get_state_count(network)[AdaptiveHypergraphs.B] == 2

    # an empty network with a certain probability of infection
    Random.seed!(42)
    network = HyperNetwork(n, 0.5, 2)
    @test get_state_count(network)[AdaptiveHypergraphs.A] == 2
    @test get_state_count(network)[AdaptiveHypergraphs.B] == 3
end

@testset "Hypernetwork: graph manipulation" begin
    n = 5
    node_state = Vector{Union{Nothing,State}}(nothing, n)
    fill!(node_state, AdaptiveHypergraphs.A)
    node_state[1] = AdaptiveHypergraphs.B
    node_state[3] = AdaptiveHypergraphs.B
    network = HyperNetwork(n, node_state, 4)

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

    # not possible to add a hyperedge of size one
    @test_throws AssertionError add_hyperedge!(network, (1,))
    # not possible to add degenerate hyperedges
    @test_throws AssertionError add_hyperedge!(network, (1, 1))
    # not possible to use nodes not in the network
    @test_throws AssertionError add_hyperedge!(network, (42, 234))
    # not possible to delete a non-existant hyperedge
    @test_throws AssertionError delete_hyperedge!(network, 213)

    include_node!(network, 3, 1)
    include_node!(network, 5, 3)
    # hyperedge 1: [1, 2, 3], hyperedge 2: [2, 3], hyperedge 3: [1, 3, 4, 5]

    # not possible to create hyperedges bigger than max_size
    @test_throws AssertionError include_node!(network, 2, 3)

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
    # get_max_size still returns the historical max value
    @test get_max_size(network) == 4

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
    network = HyperNetwork(n, 3)
    add_hyperedge!(network, (1, 3, 4))

    # test that all nodes are in the same state A (i.e., the hyperedge is not active)
    @test is_active(network, 1) == false
    @test get_num_active_hyperedges(network) == 0
    @test get_state_count(network) ==
          Dict(AdaptiveHypergraphs.A => 5, AdaptiveHypergraphs.B => 0)

    set_state!(network, 3, AdaptiveHypergraphs.B)
    # now this is no longer the case
    @test is_active(network, 1) == true
    @test get_num_active_hyperedges(network) == 1

    @test get_state_map(network) ==
          Dict(1 => AdaptiveHypergraphs.A, 2 => AdaptiveHypergraphs.A,
               3 => AdaptiveHypergraphs.B, 4 => AdaptiveHypergraphs.A,
               5 => AdaptiveHypergraphs.A)
    @test get_state_map(network, 1) ==
          Dict(1 => AdaptiveHypergraphs.A, 3 => AdaptiveHypergraphs.B,
               4 => AdaptiveHypergraphs.A)
    @test get_state_count(network) ==
          Dict(AdaptiveHypergraphs.A => 4, AdaptiveHypergraphs.B => 1)
end

@testset "Hypernetwork: graph construction" begin
    n = 50
    network = HyperNetwork(n, 0.4, 4)
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

@testset "HyperNetwork: motif count" begin
    n = 50
    network = HyperNetwork(n, 0.4, 4)
    build_RSC_hg!(network, (20, 20, 20))

    # shuffle things around
    for node in rand(1:50, 50)
        state = rand(instances(State))
        set_state!(network, node, state)
    end

    # rewire some edges to other edges
    for i in 1:60
        hyperedge = rand(get_hyperedges(network))
        node = rand(get_nodes(network, hyperedge))
        remove_node!(network, node, hyperedge)
        hyperedge_found = false
        while !hyperedge_found
            new_hyperedge = rand(get_hyperedges(network))
            if length(get_nodes(network, new_hyperedge)) == network.max_size
                continue
            end
            try
                include_node!(network, node, new_hyperedge)
            catch DegenerateHyperedge
                continue
            end
            hyperedge_found = true
        end
    end

    # rewire some edges to nodes
    for i in 1:10
        hyperedge = rand(get_hyperedges(network))
        node = rand(get_nodes(network, hyperedge))
        remove_node!(network, node, hyperedge)
        node_found = false
        while !node_found
            new_node = rand(1:n)
            if new_node == node
                continue
            end
            add_hyperedge!(network, (node, new_node))
            node_found = true
        end
    end

    # check that the number of active size-2 edges is equal to [AB]
    true_value = get_num_active_hyperedges(network, 2)
    motif_prediction = network.motif_count[Label("[AB]")]
    @test true_value == motif_prediction

    # same but for size 3
    true_value = get_num_active_hyperedges(network, 3)
    motif_prediction = network.motif_count[Label("[A2B]")] +
                       network.motif_count[Label("[AB2]")]
    @test true_value == motif_prediction

    # size 4
    true_value = get_num_active_hyperedges(network, 4)
    motif_prediction = network.motif_count[Label("[A3B]")] +
                       network.motif_count[Label("[A2B2]")] +
                       network.motif_count[Label("[AB3]")]
    @test true_value == motif_prediction

    # test that the total number of order-one motifs is equal to the number of hyperdeges
    num_order_one_motifs = 0
    for label in all_labels(4)
        if order(label) == 1
            num_order_one_motifs += network.motif_count[label]
        end
    end
    @test get_num_hyperedges(network) == num_order_one_motifs

    # count all tripples explicitly and check that we get the correct number

    # prepare the data structure
    explicit_results = Dict{Label,Int64}()
    for label in all_labels(4)
        if order(label) == 2
            explicit_results[label] = 0
        end
    end

    # COUNT THEM ALL
    for hyperedge in get_hyperedges(network)
        statecount1 = get_state_count(network, hyperedge)
        for node in get_nodes(network, hyperedge)
            int_state = get_state(network, node)
            for neighbor in get_hyperedges(network, node)
                if neighbor == hyperedge
                    continue
                end
                statecount2 = get_state_count(network, neighbor)
                label = Label(statecount1, statecount2, int_state)
                explicit_results[label] += 1
            end
        end
    end

    # the tripples are counted twice, so divide the numbers by two
    for label in keys(explicit_results)
        explicit_results[label] ÷= 2
    end

    for label in keys(explicit_results)
        println("Comparing $label with $(explicit_results[label])")
        @test explicit_results[label] == network.motif_count[label]
    end
end

# [A2B | A | A]   8 == 7    >
# [A2B | A | AB]  15 == 16  <
# [A2B | A | AB2] 11 == 12  <
# [A2B | B | B]   6 == 7    <
# [A2B | B | A2]  7 == 8    <
# [A2B | B | AB2] 7 == 8    <
# [A2B | B | A2B] 5 == 4    >
# 
# [AB2 | A | B]   8 == 9    <
# [AB2 | B | B3]  3 == 2    >
# [AB2 | B | AB2] 1 == 2    <
# [AB2 | B | AB]  7 == 6    >
# 
# [A3 | B | A]    4 == 3    >
# [A3 | A | A2]   3 == 2    >
# [A3 | A | A2B]  3 == 2    >