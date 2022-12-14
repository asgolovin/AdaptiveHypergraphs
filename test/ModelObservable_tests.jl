n = 10
node_state = Vector{Union{Nothing,State}}(nothing, n)
fill!(node_state, AdaptiveHypergraphs.S)
node_state[2] = AdaptiveHypergraphs.I
node_state[5] = AdaptiveHypergraphs.I
network = HyperNetwork(n, node_state, 4)
build_RSC_hg!(network, (3, 4, 5))

majority_voting = MajorityVoting()
rewiring_rule = RewireToSame()
adaptivity_prob = 0.5

model = DiscrModel{MajorityVoting,RewireToSame}(network,
                                                majority_voting,
                                                rewiring_rule,
                                                adaptivity_prob)

measurement_types = [StateCount, HyperedgeCount, ActiveHyperedgeCount]
mo = ModelObservable(model, measurement_types)

println(typeof.(mo.measurements))
println(measurement_types)

# check that only the required measurements were initialized
@test all(typeof.(mo.measurements) .∈ Ref(measurement_types))
@test all(measurement_types .∈ Ref(typeof.(mo.measurements)))

# check that values get recorded correctly 

record_measurements!(mo, :step)
for meas in mo.state_count
    if meas.label == S
        @test meas.values[] == [n - 2]
    else
        @test meas.values[] == [2]
    end
end

true_size_dist = [3, 4, 5]
for (i, meas) in enumerate(mo.hyperedge_count)
    @test meas.values[] == [true_size_dist[i]]
end
for (i, meas) in enumerate(mo.active_hyperedge_count)
    @test meas.values[] == [get_num_active_hyperedges(network, i + 1)]
end
