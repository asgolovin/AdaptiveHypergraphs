using AdaptiveHypergraphs

#! format: off

params = InputParams(
    NetworkParams(
        num_nodes = 100,
        num_hyperedges = (500, 500),
        infected_prob = 0.5,
    ),
    ModelParams(
        num_time_steps = Int64(1e6),
        adaptivity_rule = RewireToRandom(),
        propagation_rule = ProportionalVoting(),
        adaptivity_prob = 0.5
    ),
    VisualizationParams(
        skip_points = 1,
        buffer_size = 1000,
        misc_colormap = :Set1_7,
        panels = [:StateDistPanel,
                  :HyperedgeDistPanel,
                  :ActiveHyperedgeDistPanel,
                  :SlowManifoldPanel,
                  #:FakeDiffEqPanel,
                  :FirstOrderMotifCountPanel,
                  :SecondOrderMotifCountPanel,
                  #:MomentClosurePanel
        ]
    ),
    BatchParams(
        batch_size = 5,
        with_mpi = false,
        prompt_for_save = true,
    )
)