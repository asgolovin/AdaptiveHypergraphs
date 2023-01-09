using AdaptiveHypergraphs

#! format: off

params = InputParams(
    NetworkParams(
        num_nodes = 1000,
        num_hyperedges = (2000, 200),
        infected_prob = 0.2
    ),
    ModelParams(
        num_time_steps = Int64(2e5),
        adaptivity_rule = RewireToRandom(),
        propagation_rule = ProportionalVoting(),
        adaptivity_prob = 0.5
    ),
    VisualizationParams(
        skip_points = 100,
        buffer_size = 10000,
        panels = [:StateDistPanel,
                  :HyperedgeDistPanel,
                  :ActiveHyperedgeDistPanel,
                  :FirstOrderMotifCountPanel,
                  :SecondOrderMotifCountPanel,
                  :SlowManifoldPanel,
                  :FakeDiffEqPanel,
                  :MomentClosurePanel]
    ),
    BatchParams(
        batch_size = 4,
        with_mpi = false
    )
)