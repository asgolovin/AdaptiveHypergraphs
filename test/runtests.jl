using AdaptiveHypergraphs
using Test

@testset "AdaptiveHypergraphs.jl" begin
    include("../src/data/Network.jl")

    include("../src/simulation/AdaptivityRule.jl")
    include("../src/simulation/PropagationRule.jl")
    include("../src/simulation/Model.jl")

    include("../src/presentation/HypergraphPlot.jl")
    include("../src/presentation/ModelObservable.jl")
    include("../src/presentation/Dashboard.jl")

    @testset "HyperNetwork" begin
        @testset "HyperNetwork: constructors" begin
            # create an empty network with all suseptible nodes
            n = 5
            network = HyperNetwork(n)
            @test typeof(network.hg) <: Hypergraph
            @test size(network.hg) == (n, 0)
            @test network.hg == HyperNetwork(n).hg

            # an empty network with a given distribution of infected nodes
            node_state = Vector{Union{Nothing,State}}(nothing, n)
            fill!(node_state, S)
            node_state[2] = I
            node_state[5] = I
            network = HyperNetwork(n, node_state)
            @test network.state_dist[S] == n - 2
            @test network.state_dist[I] == 2

            # an empty network with a certain probability of infection
            Random.seed!(42)
            network = HyperNetwork(n, 0.5)
            @test network.state_dist[S] == 3
            @test network.state_dist[I] == 2
        end

        @testset "Hypernetwork: graph manipulation" begin
            n = 5
            node_state = Vector{Union{Nothing,State}}(nothing, n)
            fill!(node_state, S)
            node_state[2] = I
            node_state[5] = I
            network = HyperNetwork(n, node_state)

            add_hyperedge!(network, (1,))
            add_hyperedge!(network, (1, 2))
            add_hyperedge!(network, (2, 3))
            add_hyperedge!(network, (1, 3, 4))

            @test size(network.hg) == (n, 4)
            @test get_num_hyperedges(network) == 4
            @test get_node_degree(network, 1) == 3
            @test get_node_degree(network, 5) == 0
            hyperedge_dist = get_hyperedge_dist(network)
            @test hyperedge_dist[1] == 1
            @test hyperedge_dist[2] == 2
            @test hyperedge_dist[3] == 1

            @test_throws AssertionError add_hyperedge!(network, (42, 234))
            @test_throws AssertionError add_node!(network, (1919, 2222), S)
            @test_throws AssertionError delete_hyperedge!(network, 213)

            include_node!(network, 3, 1)
            include_node!(network, 5, 3)
            @test all([1, 3] .∈ Ref(get_nodes(network, 1)))
            @test all([2, 3, 5] .∈ Ref(get_nodes(network, 3)))
            hyperedge_dist = get_hyperedge_dist(network)
            @test hyperedge_dist[2] == 2
            @test hyperedge_dist[3] == 2

            delete_hyperedge!(network, 2)
            @test size(network.hg) == (n, 3)
            hyperedge_dist = get_hyperedge_dist(network)
            @test hyperedge_dist[2] == 1
        end

        @testset "Hypernetwork: graph info" begin
            n = 5
            network = HyperNetwork(n)
            add_hyperedge!(network, (1, 3, 4))
            # test that all nodes are in the same state S (i.e., the hyperedge is not active)
            @test is_active(network, 1) == false
            set_state!(network, 3, I)
            # now this is no longer the case
            @test is_active(network, 1) == true
        end

        @testset "Hypernetwork: graph construction" begin
            n = 50
            network = HyperNetwork(n)
            build_RSC_hg!(network, (10, 20, 30))

            @test get_num_hyperedges(network) == 60
            @test get_hyperedge_dist(network)[2] == 10
            @test get_hyperedge_dist(network)[3] == 20
            @test get_hyperedge_dist(network)[4] == 30
        end
    end

    @testset "Presentation" begin
        @testset "ModelObservable" begin
            n = 10
            node_state = Vector{Union{Nothing,State}}(nothing, n)
            fill!(node_state, S)
            node_state[2] = I
            node_state[5] = I
            network = Observable(HyperNetwork(n, node_state))
            build_RSC_hg!(network[], (3, 4, 5))

            majority_voting = MajorityVoting()
            rewiring_rule = ConflictAvoiding()
            propagation_prob = 0.5

            model = DiscrModel{MajorityVoting,ConflictAvoiding}(network[],
                                                                majority_voting,
                                                                rewiring_rule,
                                                                propagation_prob)

            mo = ModelObservable{typeof(model)}(model)
            @test typeof(mo.model) <:
                  Observable{DiscrModel{MajorityVoting,ConflictAvoiding}}
            @test mo.state_history[S][] == [n - 2]
            @test mo.state_history[I][] == [2]
            @test mo.hyperedge_history[2][] == [3]
            @test mo.hyperedge_history[3][] == [4]
            @test mo.hyperedge_history[4][] == [5]
        end
    end
end
