include("InputParams.jl")

params = InputParams(
    NetworkParams(
        num_nodes = 30,
        num_hyperedges = (30, 10, 3, 1)
    ),
    ModelParams(
        num_time_steps = 230
    ),
    VisualizationParams(
        dashboard_params = (plot_hyperedges = false, )
    )
)