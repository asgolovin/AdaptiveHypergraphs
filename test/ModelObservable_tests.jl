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
        @test series.values[] == [n - 2]
    else
        @test series.values[] == [2]
    end
end

true_size_dist = [3, 4, 5]
for (i, series) in enumerate(mo.hyperedge_series)
    @test series.values[] == [true_size_dist[i]]
end
for (i, series) in enumerate(mo.active_hyperedges_series)
    @test series.values[] == [get_num_active_hyperedges(network, i + 1)]
end