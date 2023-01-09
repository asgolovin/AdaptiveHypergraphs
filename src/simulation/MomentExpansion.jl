export moment_expansion, rhs, moment_closure

using DifferentialEquations

function moment_expansion(initial_motif_count::Dict, params, tspan)
    nparams = params.network_params

    max_size = length(nparams.num_hyperedges) + 1

    labels = filter(x -> order(x) <= 1, all_labels(max_size))
    x0 = zeros(length(labels))
    for label in labels
        x0[label_to_id(label, max_size)] = initial_motif_count[label]
    end

    problem = ODEProblem(rhs!, x0, tspan, params)

    sol = solve(problem)

    t = sol.t 
    x = sol.u

    label_to_solution = Dict{Label, Vector{Float64}}()

    solution_matrix = hcat(sol.u...)

    for label in labels
        id = label_to_id(label, max_size)
        label_to_solution[label] =  solution_matrix[id, 1:end]
    end

    return t, label_to_solution
end


function rhs!(dx, x, params, t)

    nparams = params.network_params

    max_size = length(nparams.num_hyperedges) + 1

    for i in 1:length(x)
        label = id_to_label(i, max_size)

        prop_dx = prop_term(label, x, params, max_size)
        adapt_dx = adapt_term(label, x, params, max_size)

        dx[i] = prop_dx + adapt_dx
    end

    return dx
end

function prop_term(label, x, params, max_size)
    # TODO: support for other rules. Only proportional voting is implemented at the moment! 

    # the rhs of zero-order labels is equal to zero
    if order(label) == 0
        return 0.
    end

    k = label.left_total[A]
    h = label.left_total[B]
    size = k + h

    prop_dx = 0.0
    p = params.model_params.adaptivity_prob

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
            prop_dx += n / (n + m) * (moment_closure(Label("[A$(n) B$(m-1)|B|A$(k-1)B$(h)]"), x, max_size))
            prop_dx -= m / (n + m) * (moment_closure(Label("[A$(n-1) B$(m)|A|A$(k-1)B$(h)]"), x, max_size))
        end
        if h > 0
            prop_dx -= n / (n + m) * (moment_closure(Label("[A$(n) B$(m-1)|B|A$(k)B$(h-1)]"), x, max_size))
            prop_dx += m / (n + m) * (moment_closure(Label("[A$(n-1) B$(m)|A|A$(k)B$(h-1)]"), x, max_size))
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

function adapt_term(label, x, params, max_size)
    # TODO: support for other rules. Only rewire-to-random is implemented at the moment! 

    # the rhs of zero-order labels is equal to zero
    if order(label) == 0
        return 0.
    end

    k = label.left_total[A]
    h = label.left_total[B]

    adapt_dx = 0.0
    p = params.model_params.adaptivity_prob

    # first-order term
    if k != 0 && h != 0 # active hyperedge
        adapt_dx -= x[label_to_id(label, max_size)]
    end
    
    # compute number of nodes or hyperedges which can be rewired to
    num_nodes = params.network_params.num_nodes
    small_labels = filter(x -> order(x) == 1 && size(x) != max_size, all_labels(max_size))
    num_small_hyperedges = sum([x[label_to_id(label, max_size)] for label in small_labels])
    num_candidates = num_nodes + num_small_hyperedges

    for n in 1:max_size, m in 1:(max_size - n)
        if n + m == 1
            continue
        end
        AnBm = x[label_to_id(Label("[A$(n)B$(m)]"), max_size)]

        if k > 0
            adapt_dx -= x[label_to_id(Label("[A$(k)B$(h)]"), max_size)] * AnBm / num_candidates
            adapt_dx += n / (n + m) * x[label_to_id(Label("[A$(k-1)B$(h)]"), max_size)] * AnBm / num_candidates
        end
        if h > 0
            adapt_dx += m / (n + m) * x[label_to_id(Label("[A$(k)B$(h-1)]"), max_size)] * AnBm / num_candidates
        end
    end

    if k > 0 && h > 0 && k + h < max_size
        adapt_dx += (k + 1) / (k + h + 1) * x[label_to_id(Label("[A$(k+1)B$(h)]"), max_size)]
        adapt_dx += (h + 1) / (k + h + 1) * x[label_to_id(Label("[A$(k)B$(h+1)]"), max_size)]
    end

    return adapt_dx * p
end

function moment_closure(high_order_label::Label, x::Vector{Float64}, max_size::Int64)
    @assert order(high_order_label) > 1

    left_label = Label(high_order_label.left_total)
    right_label = Label(high_order_label.right_total)
    intersection = Label(high_order_label.int)

    left_id = label_to_id(left_label, max_size)
    right_id = label_to_id(right_label, max_size)
    int_id = label_to_id(intersection, max_size)

    return x[left_id] * x[right_id] / x[int_id]
end

# the let-block immitates static variables which are saved between function calls. 
# This way, the hash table can be cached and doesn't have to be recomputed again between function calls.
let label_dict::Union{Dict, Nothing} = nothing,
    id_dict::Union{Dict, Nothing} = nothing,
    cached_max_size::Union{Int64, Nothing} = nothing

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

