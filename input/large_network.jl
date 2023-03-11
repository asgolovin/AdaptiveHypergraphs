using AdaptiveHypergraphs

#! format: off

params = InputParams(
    NetworkParams(
        num_nodes = 1000,
        num_hyperedges = (2000, 20, 20),
        state_A_prob = 0.5,
    ),
    ModelParams(
        is_discrete = true,
        num_time_steps = Int64(1e6),
        adaptivity_rule = RewireToRandom(),
        propagation_rule = ProportionalVoting(),
        adaptivity_prob = 0.5
    ),
    VisualizationParams(
        skip_points = 1000,
        buffer_size = 50000,
        misc_colormap = :Set1_7,
        panels = [:StateDistPanel,
                  :HyperedgeDistPanel,
                  :ActiveHyperedgeDistPanel,
                  :SlowManifoldPanel,
        ]
    ),
    BatchParams(
        batch_size = 5,
        with_mpi = false,
        prompt_for_save = false,
    )
)