using AdaptiveHypergraphs

include("../input/small_network.jl")

nparams, mparams, vparams = params.network_params, params.model_params, params.visualization_params

n = nparams.num_nodes
network = HyperNetwork(n, nparams.infected_prob)
build_RSC_hg!(network, nparams.num_hyperedges)

propagation_rule = mparams.propagation_rule
adaptivity_rule = mparams.adaptivity_rule

if mparams.is_discrete
    model = DiscrModel{typeof(propagation_rule), 
                       typeof(adaptivity_rule)}(network,
                                                propagation_rule,
                                                adaptivity_rule,
                                                mparams.propagation_prob)
else
    model = ContinuousModel{typeof(propagation_rule), 
                            typeof(adaptivity_rule)}(network,
                                                    propagation_rule,
                                                    adaptivity_rule,
                                                    mparams.propagation_rate,
                                                    mparams.adaptivity_rate)
end

dashboard = Dashboard(model; vparams.dashboard_params...)

if vparams.record_video
    record!(dashboard, "test_record", 100, 10, 1)
else
    run!(dashboard, mparams.num_time_steps, vparams.steps_per_update)
end