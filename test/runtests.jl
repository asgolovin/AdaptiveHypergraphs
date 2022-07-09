using AdaptiveHypergraphs
using Test

@testset "AdaptiveHypergraphs.jl" begin
    @testset "HyperNetwork" begin
        include("../src/data/Network.jl")

        @testset "HyperNetwork: constructors" begin
            # create an empty network with all suseptible nodes
            n = 5
            network = HyperNetwork(n)
            @test typeof(network.hg) <: Hypergraph
            @test size(network.hg) == (n, 0)
            @test network.hg == HyperNetwork(n).hg
            
            # an empty network with a given distribution of infected nodes
            node_state = Vector{Union{Nothing, State}}(nothing, n)
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
            node_state = Vector{Union{Nothing, State}}(nothing, n)
            fill!(node_state, S)
            node_state[2] = I
            node_state[5] = I
            network = HyperNetwork(n, node_state)

            add_hyperedge!(network, (1, ))
            add_hyperedge!(network, (1, 2))
            add_hyperedge!(network, (2, 3))
            add_hyperedge!(network, (1, 3, 4))
            
            @test size(network.hg) == (n, 4)
            @test get_num_hyperedges(network) == 4
            @test get_node_degree(network, 1) == 3
            @test get_node_degree(network, 5) == 0
        end

    end
end
