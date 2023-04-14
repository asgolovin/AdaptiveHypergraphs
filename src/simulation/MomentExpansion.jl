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

    p = mparams.adaptivity_prob
    adaptivity_rule = mparams.adaptivity_rule
    propagation_rule = mparams.propagation_rule
    num_nodes = nparams.num_nodes
    max_size = length(nparams.num_hyperedges) + 1

    if typeof(propagation_rule) <: ProportionalVoting
        convincing_prob = prop_voting_convincing_prob
    else
        convincing_prob = maj_voting_convincing_prob
    end

    if typeof(adaptivity_rule) <: RewireToRandom
        joining_prob = rtr_joining_prob
    else
        joining_prob = rts_joining_prob
    end

    for i in 1:length(x)
        motif = all_motifs(max_size)[i]
        prop_dx = prop_term(motif, x, convincing_prob, moment_closure, max_size)
        adapt_dx = adapt_term(motif, x, joining_prob, num_nodes, max_size)
        dx[i] = prop_dx * (1 - p) + adapt_dx * p
    end
    return dx
end

"""
    prop_term(motif::AbstractMotif, x::Vector, convincing_prob::Function, moment_closure::Function,
              max_size::Int64)

Compute the propagation part of the right hand-side of the ODE for motif `motif` and given
the proportional voting rule. 

# Arguments
- `motif`: the motif whos time derivative is computed
- `x`: a vector of all unknowns
- `convincing_prob`: a function `convincing_prob(state, n, m) => Float64` which gives the probabity that a node in a state `state` 
convinces the hyperedge AnBm. 
- `moment_closure`: a function which maps high-order motifs to their moment-closure. 
- `max_size`: maximum cardinality of a hyperedge
"""
function prop_term(motif::AbstractMotif, x::Vector, convincing_prob::Function,
                   moment_closure::Function,
                   max_size::Int64)
    prop_dx = 0.0

    # a set of indices of active hyperedges [A^m B^n]
    active_hyperedges = []
    for m in 1:(max_size - 1), n in 1:(max_size - m)
        push!(active_hyperedges, (m, n))
    end

    # 1. Zero-order motifs: [A] and [B]

    if order(motif) == 0
        state = motif.state
        for (m, n) in active_hyperedges
            if state == A
                prop_dx += x[index(OrderOneMotif(m, n))] *
                           (n * convincing_prob(A, m, n) - m * convincing_prob(B, m, n))
            else
                prop_dx += x[index(OrderOneMotif(m, n))] *
                           (m * convincing_prob(B, m, n) - n * convincing_prob(A, m, n))
            end
        end
        return prop_dx
    end

    # 2. First-order motifs: [A^a B^b]
    a = motif.A
    b = motif.B
    size = a + b

    # first-order terms
    if a > 0 && b > 0 # active hyperedge
        prop_dx -= x[index(motif)]
    else # inactive hyperedge
        if b == 0 # [A^a]
            # An A-node convinces [A^(a - μ) B^μ] and the hyperedge turns into [A^a]
            for mu in 1:(a - 1)
                prop_dx += convincing_prob(A, a - mu, mu) *
                           x[index(OrderOneMotif(a - mu, mu))]
            end
        elseif a == 0 # [B^b]
            # A B-node convinces [A^(nu) B^(b - nu)] and the hyperedge turns into [B^b]
            for nu in 1:(b - 1)
                prop_dx += convincing_prob(B, nu, b - nu) *
                           x[index(OrderOneMotif(nu, b - nu))]
            end
        end
    end

    # second-order terms
    for (m, n) in active_hyperedges
        # A convinces the neighboring hyperedge
        for mu in 0:m, nu in 1:n
            if mu + nu <= a
                symmetric_coeff = (m == a - nu) && (n == b + nu) ? 2.0 : 1.0
                prop_dx += symmetric_coeff * convincing_prob(A, m, n) *
                           moment_closure(OrderTwoMotif((m - mu, n - nu), (mu, nu),
                                                        (a - mu - nu, b)), x, max_size)
            end
            if mu <= a && nu <= b
                symmetric_coeff = (m == a) && (n == b) ? 2.0 : 1.0
                prop_dx -= symmetric_coeff * convincing_prob(A, m, n) *
                           moment_closure(OrderTwoMotif((m - mu, n - nu), (mu, nu),
                                                        (a - mu, b - nu)), x, max_size)
            end
        end

        # B convinces the neighboring hyperedge
        for mu in 1:m, nu in 0:n
            if mu + nu <= b
                symmetric_coeff = (m == a + mu) && (n == b - mu) ? 2.0 : 1.0
                prop_dx += symmetric_coeff * convincing_prob(B, m, n) *
                           moment_closure(OrderTwoMotif((m - mu, n - nu), (mu, nu),
                                                        (a, b - mu - nu)), x, max_size)
            end
            if mu <= a && nu <= b
                symmetric_coeff = (m == a) && (n == b) ? 2.0 : 1.0
                prop_dx -= symmetric_coeff * convincing_prob(B, m, n) *
                           moment_closure(OrderTwoMotif((m - mu, n - nu), (mu, nu),
                                                        (a - mu, b - nu)), x, max_size)
            end
        end
    end

    return prop_dx
end

function prop_voting_convincing_prob(state, m, n)
    if state == A
        return m / (m + n)
    else
        return n / (m + n)
    end
end

function maj_voting_convincing_prob(state, m, n)
    prob = state == A ? sign(m - n) : sign(n - m)
    # sign(x) returns -1, 0 or 1, so map the result to 0, 0.5 or 1
    return prob * 0.5 + 0.5
end

"""
    adapt_term(motif::AbstractMotif, x::Vector, joining_prob::Function, 
    num_nodes::Int64, max_size::Int64)

Compute the adaptivity part of the right hand-side of the ODE for motif `motif` and given
the rewire to random rule.
"""
function adapt_term(motif::AbstractMotif, x::Vector, joining_prob::Function,
                    num_nodes::Int64, max_size::Int64)
    # the rhs of zero-order motifs [A] and [B] is equal to zero
    if order(motif) == 0
        return 0.0
    end

    a = motif.A
    b = motif.B

    adapt_dx = 0.0

    # a set of indices of active hyperedges [A^m B^n]
    active_hyperedges = []
    for m in 1:(max_size - 1), n in 1:(max_size - m)
        push!(active_hyperedges, (m, n))
    end

    # Source terms: a node leaves the source hyperedge
    # The source hyperedge has to be active and has to have an allowed number of nodes

    if a + b + 1 <= max_size
        # Gain term: A leaves [A^a+1 B^b] and creates an [A^a B^b] edge
        if a + 1 > 0 && b > 0
            adapt_dx += (a + 1) / (a + b + 1) *
                        x[index(OrderOneMotif(a + 1, b))]
        end

        # Gain term: B leaves [A^a B^b+1] and creates an [A^a B^b] edge
        if a > 0 && b + 1 > 0
            adapt_dx += (b + 1) / (a + b + 1) *
                        x[index(OrderOneMotif(a, b + 1))]
        end
    end

    # Loss term: Any node leaves [A^a B^b] and destroys an [A^a B^b] edge
    if a > 0 && b > 0
        adapt_dx -= x[index(motif)]
    end

    # Target terms: a node joins the target hyperedge (or node)
    # Here, we need to sum over all possible sources A^m B^n of the node to get 
    # the probabity that, for example, an A-node was selected. 
    for (m, n) in active_hyperedges
        # The source hyperedge A^m B^n
        AmBn = x[index(OrderOneMotif(m, n))]

        # An A-node joins a hyperedge
        adapt_dx += m / (m + n) * AmBn *
                    (joining_prob(A, a - 1, b, x, max_size) -
                     joining_prob(A, a, b, x, max_size))

        # A B-node joins a hyperedge
        adapt_dx += n / (m + n) * AmBn *
                    (joining_prob(B, a, b - 1, x, max_size) -
                     joining_prob(B, a, b, x, max_size))
    end

    return adapt_dx
end

function rtr_joining_prob(state, a, b, x, max_size)
    if !(1 <= a + b < max_size && a >= 0 && b >= 0)
        return 0.0
    end

    num_nodes = x[index(OrderZeroMotif(A))] + x[index(OrderZeroMotif(B))]

    if max_size == 2
        normalization = num_nodes
    else
        # A set of motifs with less than maximum size. A new node can join those motifs. 
        small_motifs = filter(x -> order(x) == 1 && size(x) < max_size,
                              all_motifs(max_size))
        num_small_hyperedges = sum([x[index(motif)] for motif in small_motifs])
        normalization = num_nodes + num_small_hyperedges
    end

    if a == 1 && b == 0
        return x[index(OrderZeroMotif(A))] / normalization
    elseif a == 0 && b == 1
        return x[index(OrderZeroMotif(B))] / normalization
    else
        return x[index(OrderOneMotif(a, b))] / normalization
    end
end

function rts_joining_prob(state, a, b, x, max_size)
    if state == A
        if !(b == 0 && 1 <= a < max_size)
            return 0.0
        end
    else
        if !(a == 0 && 1 <= b < max_size)
            return 0.0
        end
    end

    if max_size == 2
        normalization = x[index(OrderZeroMotif(state))]
    else
        small_motifs = filter(x -> order(x) == 1 && size(x) < max_size,
                              all_motifs(max_size))
        # for rewire-to-same, the hyperedge has to additionally have all nodes in the same state
        other_state = state == A ? :B : :A
        small_inactive_motifs = filter(x -> getproperty(x, other_state) == 0,
                                       small_motifs)
        num_small_hyperedges = sum([x[index(motif)]
                                    for motif in small_inactive_motifs])
        normalization = x[index(OrderZeroMotif(state))] + num_small_hyperedges
    end

    if a == 1 && b == 0
        return x[index(OrderZeroMotif(A))] / normalization
    elseif a == 0 && b == 1
        return x[index(OrderZeroMotif(B))] / normalization
    else
        return x[index(OrderOneMotif(a, b))] / normalization
    end
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
    left_id = index(left_motif)
    right_id = index(right_motif)

    m = left_motif.A
    n = left_motif.B
    mu = high_order_motif.int.A
    nu = high_order_motif.int.B
    a = right_motif.A
    b = right_motif.B

    for (i, value) in enumerate(x)
        if isnan(value)
            println("$i, $value")
        end
    end

    numA = Int64(round(x[index(OrderZeroMotif(A))]))
    numB = Int64(round(x[index(OrderZeroMotif(B))]))

    if numA < mu || numB < nu
        return 0.0
    end

    # compute the combinatorical coefficient
    combinatorical_coeff = binomial(m, mu) * binomial(n, nu) * binomial(a, mu) *
                           binomial(b, nu)
    if left_motif == right_motif
        combinatorical_coeff /= 2
    end

    return combinatorical_coeff * x[left_id] * x[right_id] /
           (binomial(numA, mu) * binomial(numB, nu))
end