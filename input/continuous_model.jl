using AdaptiveHypergraphs

#! format: off

params = InputParams(
    NetworkParams(
        num_nodes = 1000,
        num_hyperedges = (2000, 10, 10),
        infected_prob = 0.1
    ),
    ModelParams(
        num_time_steps = 400,
        is_discrete = false,
        propagation_rate = 1,
        adaptivity_rate = 1
    ),
    VisualizationParams(
        dashboard_params = (plot_hyperedges = true, plot_hypergraph = false),
        steps_per_update = 10
    )
)