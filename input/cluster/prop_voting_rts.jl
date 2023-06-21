using AdaptiveHypergraphs

#! format: off

params = InputParams(
    NetworkParams(
        num_nodes = 10000,
        num_hyperedges = num_hyperedges = (12000, 5000, 250),
        state_A_prob = 0.5
    ),
    ModelParams(
        max_duration = 20000.0,
        adaptivity_rule = RewireToSame(),
        propagation_rule = ProportionalVoting(),    
        adaptivity_prob = vcat(collect(0.0:0.1:0.5), collect(0.525:0.025:0.7), collect(0.8:0.1:1.0)), # 17 points
    ),
    VisualizationParams(
        skip_points = 100,
        buffer_size = Int64(1e6),
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
        save_tag = "prop_voting_rts"
    )
)