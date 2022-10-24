using AdaptiveHypergraphs

#! format: off

params = InputParams(
    NetworkParams(
        num_nodes = 20,
        num_hyperedges = (20, 3, 1),
        infected_prob = 0.3
    ),
    ModelParams(
        num_time_steps = 100,
        propagation_prob = 0.2,
        adaptivity_rule = RewireToSame(),
        propagation_rule = MajorityVoting()
    ),
    VisualizationParams(
        dashboard_params = (plot_hypergraph = true,
                            plot_states = true,
                            plot_hyperedges = true,
                            plot_active_hyperedges = true),
        steps_per_update = 1
    ),
    BatchParams(
        batch_size = 10
    )
)