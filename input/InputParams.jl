using Parameters

@with_kw struct NetworkParams
    num_nodes::Integer = 100
    infected_prob::Real = 0.5
    num_hyperedges::Tuple = (100, 10)
end

@with_kw struct ModelParams
    is_discrete::Bool = true
    adaptivity_rule::AdaptivityRule = ConflictAvoiding()
    propagation_rule::PropagationRule = MajorityRule()
    num_time_steps::Integer = 500
    propagation_prob::Real = 0.5
    propagation_rate::Real = 1.
    adaptivity_rate::Real = 1.
end

@with_kw struct VisualizationParams
    record_video::Bool = false
    dashboard_params::NamedTuple = ()
    steps_per_update::Integer = 10
end

struct InputParams
    network_params::NetworkParams
    model_params::ModelParams
    visualization_params::VisualizationParams
end