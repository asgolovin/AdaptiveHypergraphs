# using AdaptiveHypergraphs
using Test

@testset "AdaptiveHypergraphs.jl" begin
    include("../src/data/Network.jl")

    include("../src/simulation/AdaptivityRule.jl")
    include("../src/simulation/PropagationRule.jl")
    include("../src/simulation/Model.jl")

    include("../src/presentation/HypergraphPlot.jl")
    include("../src/presentation/Measurements.jl")
    include("../src/presentation/ModelObservable.jl")
    include("../src/presentation/Panels.jl")
    include("../src/presentation/Dashboard.jl")
    include("../src/presentation/InputParams.jl")
    include("../src/presentation/SimulationBatch.jl")

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
            @test network.state_count[S] == n - 2
            @test network.state_count[I] == 2

            # an empty network with a certain probability of infection
            Random.seed!(42)
            network = HyperNetwork(n, 0.5)
            @test network.state_count[S] == 3
            @test network.state_count[I] == 2
        end

        @testset "Hypernetwork: graph manipulation" begin
            n = 5
            node_state = Vector{Union{Nothing,State}}(nothing, n)
            fill!(node_state, S)
            node_state[1] = I
            node_state[3] = I
            network = HyperNetwork(n, node_state)

            add_hyperedge!(network, (1, 2))
            add_hyperedge!(network, (2, 3))
            add_hyperedge!(network, (1, 3, 4))
            # hyperedge 1: [1, 2], hyperedge 2: [2, 3], hyperedge 3: [1, 3, 4]

            @test size(network.hg) == (n, 3)
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
            @test size(network.hg) == (n, 2)
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
            @test get_state_count(network) == Dict(S => 5, I => 0)

            set_state!(network, 3, I)
            # now this is no longer the case
            @test is_active(network, 1) == true
            @test get_num_active_hyperedges(network) == 1

            @test get_state_map(network) == Dict(1 => S, 2 => S, 3 => I, 4 => S, 5 => S)
            @test get_state_map(network, 1) == Dict(1 => S, 3 => I, 4 => S)
            @test get_state_count(network) == Dict(S => 4, I => 1)
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
            rewiring_rule = RewireToSame()
            propagation_prob = 0.5

            model = DiscrModel{MajorityVoting,RewireToSame}(network[],
                                                            majority_voting,
                                                            rewiring_rule,
                                                            propagation_prob)

            mo = ModelObservable{typeof(model)}(model)
            flush_buffers!(mo)
            @test typeof(mo.model) <:
                  Observable{DiscrModel{MajorityVoting,RewireToSame}}
            for series in mo.state_series
                if series.state == S
                    @test series.observable[] == [n - 2]
                else
                    @test series.observable[] == [2]
                end
            end
            for series in mo.hyperedge_series
                if series.size == 2
                    @test series.observable[] == [3]
                elseif series.size == 3
                    @test series.observable[] == [4]
                elseif series.size == 4
                    @test series.observable[] == [5]
                end
            end
            num_active = get_num_active_hyperedges(network[])
            @test mo.active_hyperedges_series.observable[] == [num_active]
        end
    end
end
