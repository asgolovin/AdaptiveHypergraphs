using AdaptiveHypergraphs

#! format: off

params = InputParams(
    NetworkParams(
        num_nodes = 1000,
        num_hyperedges = [(2000, 200), (1000, 500), (700, 700), (20, 1000)],
        state_A_prob = 0.5
    ),
    ModelParams(
        num_time_steps = Int64(5e4),
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
        with_mpi = true,
        prompt_for_save = false,
        save_tag = "saving_mpi"
    )
)