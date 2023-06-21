using AdaptiveHypergraphs

#! format: off

params = InputParams(
    NetworkParams(
        num_nodes = 1000,
        num_hyperedges = (2000, 200),
        state_A_prob = 0.5
    ),
    ModelParams(
        max_duration = 10.0,
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
                  :SlowManifoldPanel,]
    ),
    BatchParams(
        batch_size = 2,
        with_mpi = true,
        prompt_for_save = false,
        save_tag = "mpi_small_test"
    )
)