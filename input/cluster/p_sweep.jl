using AdaptiveHypergraphs

#! format: off

params = InputParams(
    NetworkParams(
        num_nodes = 10000,
        num_hyperedges = num_hyperedges = (12000, 5000, 250),
        state_A_prob = 0.5,
    ),
    ModelParams(
        max_duration = 100.0,
        adaptivity_rule = RewireToRandom(),
        propagation_rule = ProportionalVoting(),
        adaptivity_prob = vcat(collect(0.:0.2:0.6), collect(0.65:0.05:0.8), collect(0.85:0.01:1.))
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
        batch_size = 10,
        with_mpi = true,
        prompt_for_save = false,
        save_tag = "prop_voting_rtr"
    )
)