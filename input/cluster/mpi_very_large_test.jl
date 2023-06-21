using AdaptiveHypergraphs

#! format: off

params = InputParams(
    NetworkParams(
        num_nodes = Int64(1e5),
        num_hyperedges = (Int64(17e4), Int64(2e4)),
        state_A_prob = 0.5
    ),
    ModelParams(
        max_duration = 1000.0,
        adaptivity_rule = RewireToRandom(),
        propagation_rule = ProportionalVoting(),
        adaptivity_prob = 0.5
    ),
    VisualizationParams(
        skip_points = 1000,
        buffer_size = Int64(1e6),
        panels = [:StateDistPanel,
                  :HyperedgeDistPanel,
                  :ActiveHyperedgeDistPanel,
                  :FirstOrderMotifCountPanel,
                  :SlowManifoldPanel,
                  ]
    ),
    BatchParams(
        batch_size = 1,
        with_mpi = true,
        prompt_for_save = false,
        save_tag = "mpi_1e5_nodes"
    )
)