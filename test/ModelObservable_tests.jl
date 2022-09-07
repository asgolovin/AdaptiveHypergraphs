n = 10
node_state = Vector{Union{Nothing,State}}(nothing, n)
fill!(node_state, AdaptiveHypergraphs.S)
node_state[2] = AdaptiveHypergraphs.I
node_state[5] = AdaptiveHypergraphs.I
network = HyperNetwork(n, node_state)
build_RSC_hg!(network, (3, 4, 5))

majority_voting = MajorityVoting()
rewiring_rule = RewireToSame()
propagation_prob = 0.5

model = DiscrModel{MajorityVoting,RewireToSame}(network,
                                                majority_voting,
                                                rewiring_rule,
                                                propagation_prob)

mo = ModelObservable{typeof(model)}(model)
flush_buffers!(mo)
for series in mo.state_series
    if series.state == AdaptiveHypergraphs.S
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
num_active = get_num_active_hyperedges(network)
@test mo.active_hyperedges_series.observable[] == [num_active]