using AdaptiveHypergraphs

#! format: off

params = InputParams(
    NetworkParams(
        num_nodes = 1000,
        num_hyperedges = (2000, 10, 10),
        state_A_prob = 0.5
    ),
    ModelParams(
        num_time_steps = Int64(1e6),
        adaptivity_rule = RewireToRandom(),
        propagation_rule = ProportionalVoting(),
        is_discrete = false,
        propagation_rate = 1.,
        adaptivity_rate = 1.
    ),
    VisualizationParams(
        skip_points = 10,
        buffer_size = 1000
    ),
    BatchParams(
        batch_size = 10
    )
)