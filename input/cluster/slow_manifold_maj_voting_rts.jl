using AdaptiveHypergraphs

#! format: off

params = InputParams(
    NetworkParams(
        num_nodes = 10000,
        num_hyperedges = num_hyperedges = (12000, 5000, 250),
        state_A_prob = [0.2, 0.35, 0.5, 0.65, 0.8]
    ),
    ModelParams(
        max_duration = 300.0,
        adaptivity_rule = RewireToSame(),
        propagation_rule = MajorityVoting(),
        adaptivity_prob = 0.1
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
        save_tag = "slow_manifold_maj_voting_rts"
    )
)