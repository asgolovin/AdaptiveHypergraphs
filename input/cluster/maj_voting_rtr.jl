using AdaptiveHypergraphs

#! format: off

params = InputParams(
    NetworkParams(
        num_nodes = 10000,
        num_hyperedges = num_hyperedges = (12000, 5000, 250),
        state_A_prob = 0.5
    ),
    ModelParams(
        max_duration = 300.0,
        adaptivity_rule = RewireToRandom(),
        propagation_rule = MajorityVoting(),
        adaptivity_prob = collect(0.0:0.025:1.0) # 41 points
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
        save_tag = "maj_voting_rtr"
    )
)