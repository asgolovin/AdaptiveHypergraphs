export moment_expansion, rhs, moment_closure

using DifferentialEquations

function moment_expansion(params, tspan, moment_closure::Function)
    nparams = params.network_params
    num_nodes = nparams.num_nodes
    max_size = length(nparams.num_hyperedges) + 1
    A_share = nparams.infected_prob
    B_share = 1 - A_share

    # calculate the initial expected distribution of nodes and motifs
    initial_motif_count = Dict{Label,Float64}()
    for label in all_labels(max_size)
        if order(label) > 1
            continue
        end
        label_size = size(label)[1]
        if label_size == 1 # [A] or [B]
            initial_motif_count[label] = label.left[A] > 0 ? A_share * num_nodes :
                                         B_share * num_nodes
            continue
        end
        num_hyperedges = nparams.num_hyperedges[label_size - 1]
        n = label.left[A]
        m = label.left[B]
        num_AnBm = num_hyperedges * A_share^n * B_share^m * binomial(n + m, n)
        initial_motif_count[label] = num_AnBm
    end

    x0 = motif_dict_to_vector(initial_motif_count, max_size)
    p = (moment_closure, params)

    problem = ODEProblem(rhs!, x0, tspan, p)
    sol = solve(problem)
    t = sol.t

    label_to_solution = Dict{Label,Vector{Float64}}()

    solution_matrix = hcat(sol.u...)

    labels = filter(x -> order(x) <= 1, all_labels(max_size))
    for label in labels
        id = label_to_id(label, max_size)
        label_to_solution[label] = solution_matrix[id, 1:end]
    end

    return t, label_to_solution
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
        label = id_to_label(i, max_size)

        prop_dx = prop_term(propagation_rule, moment_closure, label, x, adaptivity_prob,
                            max_size)
        adapt_dx = adapt_term(adaptivity_rule, label, x, adaptivity_prob, num_nodes,
                              max_size)

        dx[i] = prop_dx + adapt_dx
    end

    return dx
end

"""
    prop_term(propagation_rule::ProportionalVoting, moment_closure, label, x,
    adaptivity_prob, max_size)

Compute the propagation part of the right hand-side of the ODE for label `label` and given
the proportional voting rule. 
"""
function prop_term(propagation_rule::ProportionalVoting, moment_closure, label, x,
                   adaptivity_prob, max_size)
    # the rhs of zero-order labels is equal to zero
    if order(label) == 0
        return 0.0
    end

    k = label.left_total[A]
    h = label.left_total[B]
    size = k + h

    prop_dx = 0.0
    p = adaptivity_prob

    # first-order term
    if k != 0 && h != 0 # active hyperedge
        prop_dx -= x[label_to_id(label, max_size)]
    else # inactive hyperedge
        for n in 1:(size - 1)
            m = size - n
            if k == 0 # [B^h]
                prop_dx += m / (n + m) * x[label_to_id(Label("[A$n B$m]"), max_size)]
            elseif h == 0 # [A^k]
                prop_dx += n / (n + m) * x[label_to_id(Label("[A$n B$m]"), max_size)]
            end
        end
    end

    # second-order terms
    for n in 1:max_size, m in 1:(max_size - n)
        if n + m == 1
            continue
        end
        if k > 0
            prop_dx += n / (n + m) *
                       (moment_closure(Label("[A$(n) B$(m-1)|B|A$(k-1)B$(h)]"), x,
                                       max_size))
            prop_dx -= m / (n + m) *
                       (moment_closure(Label("[A$(n-1) B$(m)|A|A$(k-1)B$(h)]"), x,
                                       max_size))
        end
        if h > 0
            prop_dx -= n / (n + m) *
                       (moment_closure(Label("[A$(n) B$(m-1)|B|A$(k)B$(h-1)]"), x,
                                       max_size))
            prop_dx += m / (n + m) *
                       (moment_closure(Label("[A$(n-1) B$(m)|A|A$(k)B$(h-1)]"), x,
                                       max_size))
        end
    end

    # symmetric terms
    #! format: off
    if k > 0
        prop_dx += (k - 1) / (k + h) * moment_closure(Label("[A$(k-1)B$(h)|B|A$(k-1)B$(h)]"), x, max_size)
        prop_dx -=       h / (k + h) * moment_closure(Label("[A$(k-1)B$(h)|A|A$(k-1)B$(h)]"), x, max_size)
    end
    if h > 0
        prop_dx += (h - 1) / (k + h) * moment_closure(Label("[A$(k)B$(h-1)|A|A$(k)B$(h-1)]"), x, max_size)
        prop_dx -=       k / (k + h) * moment_closure(Label("[A$(k)B$(h-1)|B|A$(k)B$(h-1)]"), x, max_size)
    end

    #! format: on
    return prop_dx * (1 - p)
end

"""
    adapt_term(adaptivity_rule::RewireToRandom, label, x, adaptivity_prob, num_nodes,
    max_size)

Compute the adaptivity part of the right hand-side of the ODE for label `label` and given
the rewire to random rule. 
"""
function adapt_term(adaptivity_rule::RewireToRandom, label, x, adaptivity_prob, num_nodes,
                    max_size)
    # TODO: support for other rules. Only rewire-to-random is implemented at the moment! 

    # the rhs of zero-order labels [A] and [B] is equal to zero
    if order(label) == 0
        return 0.0
    end

    k = label.left_total[A]
    h = label.left_total[B]

    adapt_dx = 0.0
    p = adaptivity_prob

    # Compute the number of nodes or hyperedges which can be rewired to. 
    # Those are all "small hyperedges", i.e., hyperedges with size less than the max size, 
    # plus the number of nodes. 
    if max_size == 2
        num_candidates = num_nodes
    else
        small_labels = filter(x -> order(x) == 1 && size(x)[1] < max_size,
                              all_labels(max_size))
        num_small_hyperedges = sum([x[label_to_id(label, max_size)]
                                    for label in small_labels])
        num_candidates = num_nodes + num_small_hyperedges
    end

    # Source terms: a node leaves the source hyperedge
    # The source hyperedge has to be active and has to have an allowed number of nodes

    if k + h + 1 <= max_size
        # Gain term: A leaves [A^k+1 B^h] and creates an [A^k B^h] edge
        if k + 1 > 0 && h > 0
            adapt_dx += (k + 1) / (k + h + 1) *
                        x[label_to_id(Label("[A$(k+1)B$(h)]"), max_size)]
        end

        # Gain term: B leaves [A^k B^h+1] and creates an [A^k B^h] edge
        if k > 0 && h + 1 > 0
            adapt_dx += (h + 1) / (k + h + 1) *
                        x[label_to_id(Label("[A$(k)B$(h+1)]"), max_size)]
        end
    end

    # Loss term: Any node leaves [A^k B^h] and destroys an [A^k B^h] edge
    if k > 0 && h > 0
        adapt_dx -= x[label_to_id(label, max_size)]
    end

    # Target terms: a node joins the target hyperedge (or node)
    # Here, we need to sum over all possible sources A^n B^m of the node to get 
    # the probabity that, for example, an A-node was selected. 
    for n in 1:(max_size - 1), m in 1:(max_size - n)
        # The source hyperedge A^n B^m
        AnBm = x[label_to_id(Label("[A$(n)B$(m)]"), max_size)]

        # Gain term: an A-node joins A^k-1 B^h and creates an A^k B^h edge
        if k > 0
            adapt_dx += n / (n + m) *
                        x[label_to_id(Label("[A$(k-1)B$(h)]"), max_size)] *
                        AnBm / num_candidates
        end

        # Gain term: a B-node joins A^k B^h-1 and creates an A^k B^h edge. 
        if h > 0
            adapt_dx += m / (n + m) *
                        x[label_to_id(Label("[A$(k)B$(h-1)]"), max_size)] *
                        AnBm / num_candidates
        end

        # Loss term: any node joins A^k B^h and destroys the A^k B^h edge
        if k + h < max_size
            adapt_dx -= x[label_to_id(Label("[A$(k)B$(h)]"), max_size)] * AnBm /
                        num_candidates
        end
    end

    return adapt_dx * p
end

function adapt_term(adaptivity_rule::RewireToSame, label, x, adaptivity_prob, num_nodes,
                    max_size)
    # NOTE: NOT TESTED, CAN BE WRONG

    # the rhs of zero-order labels [A] and [B] is equal to zero
    if order(label) == 0
        return 0.0
    end

    k = label.left_total[A]
    h = label.left_total[B]

    adapt_dx = 0.0
    p = adaptivity_prob

    num_candidates = Dict{State,Float64}()
    # Compute the number of nodes or hyperedges which can be rewired to. 
    # Those are all "small hyperedges", i.e., hyperedges with size less than the max size, 
    # plus the number of nodes. 
    if max_size == 2
        num_candidates[A] = x[label_to_id(Label("[A]"), max_size)]
        num_candidates[B] = x[label_to_id(Label("[B]"), max_size)]
    else
        # TODO TODO TODO CHANGE CHANGE CHANGE
        small_labels = filter(x -> order(x) == 1 && size(x)[1] < max_size,
                              all_labels(max_size))
        num_small_hyperedges = sum([x[label_to_id(label, max_size)]
                                    for label in small_labels])
        num_candidates = num_nodes + num_small_hyperedges
    end

    # Source terms: a node leaves the source hyperedge
    # Only possible for active hyperedges, i.e., for k, h > 0.
    if k > 0 && h > 0
        if k + h < max_size
            # Gain term: A leaves [A^k+1 B^h] and creates an [A^k B^h] edge
            adapt_dx += (k + 1) / (k + h + 1) *
                        x[label_to_id(Label("[A$(k+1)B$(h)]"), max_size)]
            # Gain term: B leaves [A^k B^h+1] and creates an [A^k B^h] edge
            adapt_dx += (h + 1) / (k + h + 1) *
                        x[label_to_id(Label("[A$(k)B$(h+1)]"), max_size)]
        end
        # Loss term: Any node leaves [A^k B^h] and destroys an [A^k B^h] edge
        adapt_dx -= x[label_to_id(label, max_size)]
    end

    # Target terms: a node joins the target hyperedge (or node)
    # Here, we need to sum over all possible sources A^n B^m of the node to get 
    # the probabity that, for example, an A-node was selected. 
    # For rewire-to-same, a node may only join a hyperedge if *all* nodes in this 
    # hyperedge have the same opinion. 
    for n in 1:(max_size - 1), m in 1:(max_size - n)
        # The source hyperedge A^n B^m
        AnBm = x[label_to_id(Label("[A$(n)B$(m)]"), max_size)]

        # Gain term: an A-node joins A^k-1 B^h and creates an A^k B^h edge
        if k > 0 && h == 0
            adapt_dx += n / (n + m) *
                        x[label_to_id(Label("[A$(k-1)B$(h)]"), max_size)] *
                        AnBm / num_candidates[A]
        end

        # Gain term: a B-node joins A^k B^h-1 and creates an A^k B^h edge. 
        if h > 0 && k == 0
            adapt_dx += m / (n + m) *
                        x[label_to_id(Label("[A$(k)B$(h-1)]"), max_size)] *
                        AnBm / num_candidates[B]
        end

        # Loss term: an A-node joins A^k B^h and destroys the A^k B^h edge
        if k + h < max_size && h == 0
            adapt_dx -= x[label_to_id(Label("[A$(k)B$(h)]"), max_size)] * AnBm /
                        num_candidates
        end
        # Loss term: an A-node joins A^k B^h and destroys the A^k B^h edge
        if k + h < max_size && h == 0
            adapt_dx -= n / (n + m) * x[label_to_id(Label("[A$(k)B$(h)]"), max_size)] *
                        AnBm / num_candidates[B]
        end
        # Loss term: a B-node joins A^k B^h and destroys the A^k B^h edge
        if k + h < max_size && k == 0
            adapt_dx -= m / (n + m) * x[label_to_id(Label("[A$(k)B$(h)]"), max_size)] *
                        AnBm / num_candidates[B]
        end
    end

    return adapt_dx * p
end

function motif_dict_to_vector(motif_count::Dict, max_size::Int64)
    labels = filter(x -> order(x) <= 1, all_labels(max_size))
    x = zeros(length(labels))
    for label in labels
        x[label_to_id(label, max_size)] = motif_count[label]
    end
    return x
end

function moment_closure(high_order_label::Label, x::Vector{Float64}, max_size::Int64)
    @assert order(high_order_label) > 1

    left_label = Label(high_order_label.left_total)
    right_label = Label(high_order_label.right_total)
    intersection = Label(high_order_label.int)

    left_id = label_to_id(left_label, max_size)
    right_id = label_to_id(right_label, max_size)
    int_id = label_to_id(intersection, max_size)

    # compute the combinatorical prefactor
    int_state = high_order_label.int_state
    left_count = high_order_label.left_total[int_state]
    right_count = high_order_label.right_total[int_state]
    issymmetrical = left_label == right_label ? 0.5 : 1
    prefactor = issymmetrical * left_count * right_count

    return prefactor * x[left_id] * x[right_id] / x[int_id]
end

# the let-block immitates static variables which are saved between function calls. 
# This way, the hash table can be cached and doesn't have to be recomputed again between function calls.
let label_dict::Union{Dict,Nothing} = nothing,
    id_dict::Union{Dict,Nothing} = nothing,
    cached_max_size::Union{Int64,Nothing} = nothing

    """
    Return the index of the label in the vector x if the maximum size of a hyperedge is `max_size`.
    """
    global function label_to_id(label::Label, max_size::Int64)
        # compute the hash table for the first time
        if isnothing(label_dict) || cached_max_size != max_size
            labels = filter(x -> order(x) <= 1, all_labels(max_size))
            label_dict = Dict(label => i for (i, label) in enumerate(labels))
            id_dict = Dict(i => label for (i, label) in enumerate(labels))
            cached_max_size = max_size
        end
        return label_dict[label]
    end

    """
    Return the label corresponding to an index in the vector x if the maximum size of a hyperedge is `max_size`.
    """
    global function id_to_label(id::Int64, max_size::Int64)
        if isnothing(id_dict) || cached_max_size != max_size
            label_to_id(Label("[ A ]"), max_size)
        end
        return id_dict[id]
    end
end
