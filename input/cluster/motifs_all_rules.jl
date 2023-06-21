using AdaptiveHypergraphs

#! format: off

params = InputParams(
    NetworkParams(
        num_nodes = 10000,
        num_hyperedges = (12000, 5000, 250),
        state_A_prob = 0.75,
    ),
    ModelParams(
        max_duration = 100.0,
        adaptivity_rule = [RewireToRandom(), RewireToSame()],
        propagation_rule = [ProportionalVoting(), MajorityVoting()],
        adaptivity_prob = 0.5
    ),
    VisualizationParams(
        skip_points = 1000,
        buffer_size = Int64(1e5),
        panels = [:StateDistPanel,
                  :HyperedgeDistPanel,
                  :ActiveHyperedgeDistPanel,
                  :SlowManifoldPanel,
                  :FirstOrderMotifCountPanel
                  ]
    ),
    BatchParams(
        batch_size = 12,
        with_mpi = true,
        prompt_for_save = false,
        save_tag = "motifs_all_rules"
    )
)