using AdaptiveHypergraphs

include("../input/small_network.jl")

nparams, mparams, vparams = params.network_params, params.model_params, params.visualization_params

n = nparams.num_nodes
network = HyperNetwork(n, 0.5)
build_RSC_hg!(network, nparams.num_hyperedges)

propagation_rule = mparams.propagation_rule
adaptivity_rule = mparams.adaptivity_rule

if mparams.is_discrete
    model_type = DiscrModel{typeof(propagation_rule), typeof(adaptivity_rule)}
else
    model_type = ContModel{typeof(propagation_rule), typeof(adaptivity_rule)}
end

model = model_type(network, propagation_rule, adaptivity_rule)
@show vparams.dashboard_params
dashboard = Dashboard(model; vparams.dashboard_params...)

if !vparams.record_video
    run!(dashboard, mparams.num_time_steps, vparams.steps_per_update)
else
    record!(dashboard, "test_record", 100, 10, 1)
end