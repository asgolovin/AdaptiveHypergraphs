using AdaptiveHypergraphs

intersection_size = Dict(1 => 0, 2 => 0, 3 => 0)
n = 500
size = (1000, 1000, 1000)

for b in 1:100
    local network = HyperNetwork(n, 0.5)
    build_RSC_hg!(network, size)

    for he in 1:get_num_hyperedges(network)
        first_nodes = get_nodes(network, he)
        for node in first_nodes
            neighbors = get_hyperedges(network, node)
            filter!(x -> x != he, neighbors)
            for neighbor in neighbors
                second_nodes = get_nodes(network, neighbor)
                int_size = length(intersect(first_nodes, second_nodes))
                intersection_size[int_size] += 1
            end
        end
    end
end

share = intersection_size[2] / intersection_size[1]
println("share of 2-intersections vs 1-intersection: $(share * 100)%")