using AdaptiveHypergraphs

#! format: off

params = InputParams(
    NetworkParams(
        num_nodes = 20,
        num_hyperedges = (20, ),
        state_A_prob = 0.3
    ),
    ModelParams(
        max_duration = 100.0,
        adaptivity_prob = 0.2,
        adaptivity_rule = RewireToSame(),
        propagation_rule = MajorityVoting()
    ),
    VisualizationParams(
        skip_points = 3,
        buffer_size = 5
    ),
    BatchParams(
        batch_size = 10
    )
)