using AdaptiveHypergraphs

#! format: off

params = InputParams(
    NetworkParams(
        num_nodes = 1000,
        num_hyperedges = (200, 1200, 20),
        infected_prob = 0.5
    ),
    ModelParams(
        num_time_steps = Int64(1e8),
        adaptivity_rule = RewireToRandom(),
        propagation_rule = ProportionalVoting(),
        adaptivity_prob = collect(0.7:0.05:1.)
    ),
    VisualizationParams(
        skip_points = 100,
        buffer_size = 100000
    ),
    BatchParams(
        batch_size = 1
    )
)