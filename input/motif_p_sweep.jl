using AdaptiveHypergraphs

#! format: off

params = InputParams(
    NetworkParams(
        num_nodes = 10000,
        num_hyperedges = (12000, 5000, 250),
        state_A_prob = 0.5,
    ),
    ModelParams(
        is_discrete = true,
        max_duration = 10.0,
        adaptivity_rule = RewireToRandom(),
        propagation_rule = MajorityVoting(),
        adaptivity_prob = collect(0.0:0.1:1.0)
    ),
    VisualizationParams(
        skip_points = 100,
        buffer_size = Int64(1e5),
        panels = [:StateDistPanel,
            :HyperedgeDistPanel,
            :ActiveHyperedgeDistPanel,
            :SlowManifoldPanel,
            :FirstOrderMotifCountPanel,
            :SecondOrderMotifCountPanel,
        ]
    ),
    BatchParams(
        batch_size = 12,
        with_mpi = true,
        prompt_for_save = false,
        save_tag = "motifs_p_sweep_maj_rtr",
    )
)