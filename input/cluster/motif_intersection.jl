using AdaptiveHypergraphs

#! format: off

params = InputParams(
    NetworkParams(
        num_nodes = 1000,
        num_hyperedges = num_hyperedges = (5000, 5000, 5000),
        state_A_prob = 0.5,
    ),
    ModelParams(
        max_duration = 100.0,
        adaptivity_rule = RewireToRandom(),
        propagation_rule = ProportionalVoting(),
        adaptivity_prob = 0.5
    ),
    VisualizationParams(
        skip_points = 10,
        buffer_size = Int64(1e3),
        panels = [:StateDistPanel,
                  :HyperedgeDistPanel,
                  :ActiveHyperedgeDistPanel,
                  :SlowManifoldPanel,
                  :FirstOrderMotifCountPanel,
                  :SecondOrderMotifCountPanel,
                  ]
    ),
    BatchParams(
        batch_size = 5,
        with_mpi = true,
        prompt_for_save = false,
        save_tag = "motif_intersection_N_1e3"
    )
)