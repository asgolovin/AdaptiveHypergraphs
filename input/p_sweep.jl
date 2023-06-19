using AdaptiveHypergraphs

#! format: off

params = InputParams(
    NetworkParams(
        num_nodes = 1000,
        num_hyperedges = (1200, 500, 25),
        state_A_prob = 0.5,
    ),
    ModelParams(
        is_discrete = true,
        max_duration = 100.,
        adaptivity_rule = [RewireToSame(), RewireToRandom()],
        propagation_rule = [MajorityVoting(), ProportionalVoting()],
        adaptivity_prob = collect(0.0:0.05:1.0),
    ),
    VisualizationParams(
        skip_points = 100,
        buffer_size = 10000,
        misc_colormap = :Set1_7,
        panels = [:StateDistPanel,
                  :HyperedgeDistPanel,
                  :ActiveHyperedgeDistPanel,
                  :SlowManifoldPanel,
        ]
    ),
    BatchParams(
        batch_size = 5,
        with_mpi = true,
        prompt_for_save = true,
        save_tag = "p_sweep_small_network",
    )
)