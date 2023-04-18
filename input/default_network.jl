using AdaptiveHypergraphs

#! format: off

params = InputParams(
    NetworkParams(
        num_nodes = 1000,
        num_hyperedges = (2000, 20, 20),
        state_A_prob = 0.75,
    ),
    ModelParams(
        is_discrete = true,
        max_duration = 100.,
        adaptivity_rule = RewireToSame(),
        propagation_rule = ProportionalVoting(),
        adaptivity_prob = 0.55
    ),
    VisualizationParams(
        skip_points = 1000,
        buffer_size = 10000,
        misc_colormap = :Set1_7,
        panels = [:StateDistPanel,
                  :HyperedgeDistPanel,
                  :ActiveHyperedgeDistPanel,
                  :SlowManifoldPanel,
        ]
    ),
    BatchParams(
        batch_size = 3,
        with_mpi = false,
        prompt_for_save = false,
    )
)