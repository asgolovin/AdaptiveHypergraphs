using AdaptiveHypergraphs

#! format: off

params = InputParams(
    NetworkParams(
        num_nodes = 10000,
        num_hyperedges = num_hyperedges = (12000, 5000, 250),
        state_A_prob = [0.2, 0.35, 0.5, 0.65, 0.8]
    ),
    ModelParams(
        max_duration = 5000.0,
        adaptivity_rule = RewireToRandom(),
        propagation_rule = ProportionalVoting(),
        adaptivity_prob = 0.7
    ),
    VisualizationParams(
        skip_points = 100,
        buffer_size = Int64(1e5),
        panels = [:StateDistPanel,
                  :HyperedgeDistPanel,
                  :ActiveHyperedgeDistPanel,
                  :SlowManifoldPanel,
                  ]
    ),
    BatchParams(
        batch_size = 12,
        with_mpi = true,
        prompt_for_save = false,
        save_tag = "slow_manifold_prop_voting_rtr"
    )
)