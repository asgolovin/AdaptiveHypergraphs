using AdaptiveHypergraphs


n = 100
network = HyperNetwork(n, 0.5)
build_RSC_hg!(network, (2n, n ÷ 2, 10, 1))

majority_rule = MajorityRule()
rewiring_rule = RewiringRule(0.5)

model = DiscrModel{MajorityRule, RewiringRule}(network,
                                               majority_rule,
                                               rewiring_rule)

dashboard = Dashboard(model; plot_hypergraph=true, is_interactive=false)
run!(dashboard, 500, 10)
# record!(dashboard, "test_record", 100, 10, 1)