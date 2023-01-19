using AdaptiveHypergraphs

#! format: off

params = InputParams(
    NetworkParams(
        num_nodes = 1000,
        num_hyperedges = (2000, 200),
        infected_prob = 0.5
    ),
    ModelParams(
        num_time_steps = Int64(2e5),
        adaptivity_rule = RewireToRandom(),
        propagation_rule = ProportionalVoting(),
        adaptivity_prob = 0.5
    ),
    VisualizationParams(
        skip_points = 500,
        buffer_size = 5000,
        misc_colormap = :Set1_7,
        panels = [#:StateDistPanel,
                  #:HyperedgeDistPanel,
                  #:ActiveHyperedgeDistPanel,
                  #:SlowManifoldPanel,
                  :FakeDiffEqPanel,
                  #:SecondOrderMotifCountPanel,
                  :MomentClosurePanel
        ]
    ),
    BatchParams(
        batch_size = 10,
        with_mpi = false
    )
)