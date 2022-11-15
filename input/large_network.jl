using AdaptiveHypergraphs

#! format: off

params = InputParams(
    NetworkParams(
        num_nodes = 1000,
        num_hyperedges = (2000, 200, 20)
    ),
    ModelParams(
        num_time_steps = Int64(1e6),
        adaptivity_rule = RewireToRandom(),
        propagation_rule = ProportionalVoting(),
        adaptivity_prob = 0.5
    ),
    VisualizationParams(
        skip_points = 100,
        buffer_size = 100000
    ),
    BatchParams(
        batch_size = 10
    )
)