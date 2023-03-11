export moment_expansion, rhs, moment_closure

using DifferentialEquations

"""
    moment_expansion(params::InputParams, tspan::NTuple{2,Float64}, moment_closure::Function)

Compute the mean-field time evolution of the system specified by `params`.

# Arguments
- `params`: the InputParams used to specify the model and the network
- `tspan`: the time interval on which to simulate the system
- `moment_closure`: a funciton which expresses high-order motifs in terms of low-order ones

# Return
- `t`: a vector of time points
- `motif_to_solution`: a dict which maps first-order motifs to the time evolution of the corresponding motif
"""
function moment_expansion(params::InputParams, tspan::NTuple{2,Float64},
                          moment_closure::Function)
    nparams = params.network_params
    num_nodes = nparams.num_nodes
    max_size = length(nparams.num_hyperedges) + 1
    A_share = nparams.state_A_prob
    B_share = 1 - A_share

    # calculate the initial expected distribution of nodes and motifs
    initial_motif_count = Dict{AbstractMotif,Float64}()
    for motif in all_motifs(max_size)
        if order(motif) == 0 # [A] or [B]
            initial_motif_count[motif] = motif.state == A ? A_share * num_nodes :
                                         B_share * num_nodes
            continue
        elseif order(motif) == 2
            continue
        end

        motif_size = size(motif)
        num_hyperedges = nparams.num_hyperedges[motif_size - 1]
        n = motif.A
        m = motif.B
        num_AnBm = num_hyperedges * A_share^n * B_share^m * binomial(n + m, n)
        initial_motif_count[motif] = num_AnBm
    end

    x0 = motif_dict_to_vector(initial_motif_count, max_size)
    p = (moment_closure, params)

    problem = ODEProblem(rhs!, x0, tspan, p)
    sol = solve(problem)
    t = sol.t

    motif_to_solution = Dict{AbstractMotif,Vector{Float64}}()

    solution_matrix = hcat(sol.u...)

    motifs = filter(x -> order(x) <= 1, all_motifs(max_size))
    for motif in motifs
        id = index(motif)
        motif_to_solution[motif] = solution_matrix[id, 1:end]
    end

    return t, motif_to_solution
end

function rhs!(dx, x, p, t)
    moment_closure, params = p

    nparams = params.network_params
    mparams = params.model_params

    adaptivity_prob = mparams.adaptivity_prob
    adaptivity_rule = mparams.adaptivity_rule
    propagation_rule = mparams.propagation_rule
    num_nodes = nparams.num_nodes
    max_size = length(nparams.num_hyperedges) + 1

    for i in 1:length(x)
        motif = all_motifs(max_size)[i]

        prop_dx = prop_term(propagation_rule, moment_closure, motif, x, adaptivity_prob,
                            max_size)
        adapt_dx = adapt_term(adaptivity_rule, motif, x, adaptivity_prob, num_nodes,
                              max_size)

        dx[i] = prop_dx + adapt_dx
    end

    return dx
end

"""
    prop_term(propagation_rule::ProportionalVoting, moment_closure, motif, x,
    adaptivity_prob, max_size)

Compute the propagation part of the right hand-side of the ODE for motif `motif` and given
the proportional voting rule. 
"""
function prop_term(propagation_rule::ProportionalVoting, moment_closure, motif, x,
                   adaptivity_prob, max_size)
    # the rhs of zero-order motifs is equal to zero
    if order(motif) == 0
        return 0.0
    end

    k = motif.A
    h = motif.B
    size = k + h

    prop_dx = 0.0
    p = adaptivity_prob

    # first-order term
    if k != 0 && h != 0 # active hyperedge
        prop_dx -= x[index(motif)]
    else # inactive hyperedge
        for n in 1:(size - 1)
            m = size - n
            if k == 0 # [B^h]
                prop_dx += m / (n + m) * x[index(OrderOneMotif(n, m))]
            elseif h == 0 # [A^k]
                prop_dx += n / (n + m) * x[index(OrderOneMotif(n, m))]
            end
        end
    end

    # second-order terms
    for n in 1:max_size, m in 1:(max_size - n)
        if k > 0
            prop_dx += n / (n + m) *
                       (moment_closure(OrderTwoMotif((n, m - 1), (0, 1), (k - 1, h)), x,
                                       max_size))
            prop_dx -= m / (n + m) *
                       (moment_closure(OrderTwoMotif((n - 1, m), (1, 0), (k - 1, h)), x,
                                       max_size))
        end
        if h > 0
            prop_dx -= n / (n + m) *
                       (moment_closure(OrderTwoMotif((n, m - 1), (0, 1), (k, h - 1)), x,
                                       max_size))
            prop_dx += m / (n + m) *
                       (moment_closure(OrderTwoMotif((n - 1, m), (1, 0), (k, h - 1)), x,
                                       max_size))
        end
    end

    # symmetric terms
    #! format: off
    if k > 0
        prop_dx += (k - 1) / (k + h) * moment_closure(OrderTwoMotif((k-1, h), (0, 1), (k-1, h)), x, max_size)
        prop_dx -=       h / (k + h) * moment_closure(OrderTwoMotif((k-1, h), (1, 0), (k-1, h)), x, max_size)
    end
    if h > 0
        prop_dx += (h - 1) / (k + h) * moment_closure(OrderTwoMotif((k, h-1), (1, 0), (k, h-1)), x, max_size)
        prop_dx -=       k / (k + h) * moment_closure(OrderTwoMotif((k, h-1), (0, 1), (k, h-1)), x, max_size)
    end

    #! format: on
    return prop_dx * (1 - p)
end

"""
    adapt_term(adaptivity_rule::RewireToRandom, motif, x, adaptivity_prob, num_nodes,
    max_size)

Compute the adaptivity part of the right hand-side of the ODE for motif `motif` and given
the rewire to random rule. 
"""
function adapt_term(adaptivity_rule::RewireToRandom, motif, x, adaptivity_prob,
                    num_nodes, max_size)
    # TODO: support for other rules. Only rewire-to-random is implemented at the moment! 

    # the rhs of zero-order motifs [A] and [B] is equal to zero
    if order(motif) == 0
        return 0.0
    end

    k = motif.A
    h = motif.B

    adapt_dx = 0.0
    p = adaptivity_prob

    # Compute the number of nodes or hyperedges which can be rewired to. 
    # Those are all "small hyperedges", i.e., hyperedges with size less than the max size, 
    # plus the number of nodes. 
    if max_size == 2
        num_candidates = num_nodes
    else
        small_motifs = filter(x -> order(x) == 1 && size(x)[1] < max_size,
                              all_motifs(max_size))
        num_small_hyperedges = sum([x[index(motif)] for motif in small_motifs])
        num_candidates = num_nodes + num_small_hyperedges
    end

    # Source terms: a node leaves the source hyperedge
    # The source hyperedge has to be active and has to have an allowed number of nodes

    if k + h + 1 <= max_size
        # Gain term: A leaves [A^k+1 B^h] and creates an [A^k B^h] edge
        if k + 1 > 0 && h > 0
            adapt_dx += (k + 1) / (k + h + 1) *
                        x[index(OrderOneMotif(k + 1, h))]
        end

        # Gain term: B leaves [A^k B^h+1] and creates an [A^k B^h] edge
        if k > 0 && h + 1 > 0
            adapt_dx += (h + 1) / (k + h + 1) *
                        x[index(OrderOneMotif(k, h + 1))]
        end
    end

    # Loss term: Any node leaves [A^k B^h] and destroys an [A^k B^h] edge
    if k > 0 && h > 0
        adapt_dx -= x[index(motif)]
    end

    # Target terms: a node joins the target hyperedge (or node)
    # Here, we need to sum over all possible sources A^n B^m of the node to get 
    # the probabity that, for example, an A-node was selected. 
    for n in 1:(max_size - 1), m in 1:(max_size - n)
        # The source hyperedge A^n B^m
        AnBm = x[index(OrderOneMotif(n, m))]

        # Gain term: an A-node joins A^k-1 B^h and creates an A^k B^h edge
        if k == 1 && h == 1
            adapt_dx += n / (n + m) *
                        x[index(OrderZeroMotif(B))] * AnBm / num_candidates
        elseif k == 2 && h == 0
            adapt_dx += n / (n + m) *
                        x[index(OrderZeroMotif(A))] * AnBm / num_candidates
        elseif k > 0
            adapt_dx += n / (n + m) *
                        x[index(OrderOneMotif(k - 1, h))] * AnBm / num_candidates
        end

        # Gain term: a B-node joins A^k B^h-1 and creates an A^k B^h edge. 
        if k == 1 && h == 1
            adapt_dx += m / (n + m) *
                        x[index(OrderZeroMotif(A))] * AnBm / num_candidates
        elseif k == 0 && h == 2
            adapt_dx += m / (n + m) *
                        x[index(OrderZeroMotif(B))] * AnBm / num_candidates
        elseif h > 0
            adapt_dx += m / (n + m) *
                        x[index(OrderOneMotif(k, h - 1))] * AnBm / num_candidates
        end

        # Loss term: any node joins A^k B^h and destroys the A^k B^h edge
        if k + h < max_size
            adapt_dx -= x[index(OrderOneMotif(k, h))] * AnBm / num_candidates
        end
    end

    return adapt_dx * p
end

function motif_dict_to_vector(motif_count::Dict, max_size::Int64)
    motifs = filter(x -> order(x) <= 1, all_motifs(max_size))
    x = zeros(length(motifs))
    for motif in motifs
        x[index(motif)] = motif_count[motif]
    end
    return x
end

function moment_closure(high_order_motif::OrderTwoMotif, x::Vector{Float64},
                        max_size::Int64)
    @assert order(high_order_motif) > 1

    left_motif = high_order_motif.left_motif
    right_motif = high_order_motif.right_motif
    int_state = high_order_motif.int.A > 0 ? A : B

    left_id = index(left_motif)
    right_id = index(right_motif)
    int_id = index(OrderZeroMotif(int_state))

    # compute the combinatorical prefactor
    if int_state == A
        left_count = high_order_motif.left_motif.A
        right_count = high_order_motif.right_motif.A
    else
        left_count = high_order_motif.left_motif.B
        right_count = high_order_motif.right_motif.B
    end
    issymmetrical = left_motif == right_motif ? 0.5 : 1
    prefactor = issymmetrical * left_count * right_count

    return prefactor * x[left_id] * x[right_id] / x[int_id]
end