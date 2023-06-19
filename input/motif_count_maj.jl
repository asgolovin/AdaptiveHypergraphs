using AdaptiveHypergraphs

#! format: off

params = InputParams(
    NetworkParams(
        num_nodes = 10000,
        num_hyperedges = (12000, 5000, 250),
        state_A_prob = 0.75,
    ),
    ModelParams(
        is_discrete = true,
        max_duration = 100.,
        adaptivity_rule = [RewireToSame(), RewireToRandom()],
        propagation_rule = MajorityVoting(),
        adaptivity_prob = 0.1,
    ),
    VisualizationParams(
        skip_points = 100,
        buffer_size = 10000,
        misc_colormap = :Set1_7,
        panels = [:StateDistPanel,
                  :HyperedgeDistPanel,
                  :ActiveHyperedgeDistPanel,
                  :SlowManifoldPanel,
                  :FirstOrderMotifCountPanel,
        ]
    ),
    BatchParams(
        batch_size = 12,
        with_mpi = true,
        prompt_for_save = false,
        save_tag = "motifs_maj",
    )
)