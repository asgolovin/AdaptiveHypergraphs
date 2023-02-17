using AdaptiveHypergraphs

#! format: off

params = InputParams(
    NetworkParams(
        num_nodes = 20,
        num_hyperedges = (20, ),
        infected_prob = 0.3
    ),
    ModelParams(
        num_time_steps = 100,
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