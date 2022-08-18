using AdaptiveHypergraphs

include("../input/small_network.jl")

function _create_model(network, mparams)
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
    return model
end

nparams, mparams, vparams, bparams = params.network_params, 
                                     params.model_params,
                                     params.visualization_params,
                                     params.batch_params

n = nparams.num_nodes
network = HyperNetwork(n, nparams.infected_prob)
build_RSC_hg!(network, nparams.num_hyperedges)

model = _create_model(network, mparams)

dashboard = Dashboard(model; vparams.dashboard_params...)

for t in 1:bparams.batch_size
    reset!(dashboard, model)
    run!(dashboard, mparams.num_time_steps, vparams.steps_per_update)
    global network = HyperNetwork(n, nparams.infected_prob)
    build_RSC_hg!(network, nparams.num_hyperedges)
    global model = _create_model(network, mparams)
end