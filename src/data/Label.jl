using AutoHashEquals

export Hyperedge, Label, State, order, size, all_labels

@enum State::Bool A = false B = true

const Hyperedge = Dict{State,Int64}

@auto_hash_equals struct Label
    left::Hyperedge
    int::Hyperedge
    right::Hyperedge

    function Label(left::Hyperedge, int::Hyperedge, right::Hyperedge)
        # check that all dicts have all keys
        for state in instances(State), hyperedge in [left, int, right]
            if state ∉ keys(hyperedge)
                hyperedge[state] = 0
            end
        end

        # Switch the order of the left and right hyperedges such that 
        # the hyperedge with the larger size comes first. 
        # If both hyperedges have equal size, choose the one with more A-nodes. 
        size_left = left[A] + left[B]
        size_right = right[A] + right[B]

        if size_left > size_right
            return new(left, int, right)
        elseif size_right > size_left
            return new(right, int, left)
        end

        if left[A] >= right[A]
            return new(left, int, right)
        end

        return new(right, int, left)
    end
end

"""
    Label(str::String)

Construct labels from string notation [AB|A|A3B].

Possible input options:
 - [`<patt>`] for a single hyperedge
 - [`<patt1>`|`<patt2>`|`<patt3>`] for an intersection of two hyperedges

The patterns can be in any form of:
 - `A3 B2` - three A-nodes and two B-nodes.
 - `A B2` - one A-node and two B-nodes. Missing number is interpreted as one node.
 - `B2` - zero A-nodes and two B-nodes. Missing letter is interpreted as zero nodes.

# Examples:
    julia> Label("[AB]")
    Label(Dict(A => 1, B => 1), Dict(A => 0, B => 0), Dict(A => 0, B => 0))

    julia> Label("[AB2 | B | A2B3]")
    Label(Dict(A => 1, B => 2), Dict(A => 0, B => 1), Dict(A => 2, B => 3))

"""
function Label(str::String)
    inner_patt = r"(A(\d*))?(B(\d*))?"
    outter_patt = r"\[([^\|]*)(\|(.*)\|(.*))?\]"

    # remove whitespaces
    str = replace(str, r"\s" => "")

    outter_match = match(outter_patt, str)

    if outter_match === nothing
        throw(ArgumentError("The label $str could not be parsed: the outter pattern could not be matched"))
    end

    inner_patterns = [outter_match[1], outter_match[3], outter_match[4]]

    hyperedges = []

    for i in 1:3
        if inner_patterns[i] === nothing
            if i == 1
                throw(ArgumentError("The label $str could not be parsed: the first pattern not found."))
            else
                hyperedge = Dict(A => 0, B => 0)
                push!(hyperedges, hyperedge)
                continue
            end
        end

        inner_match = match(inner_patt, inner_patterns[i])
        if inner_match === nothing
            throw(ArgumentError("The label $str could not be parsed: subpattern $(inner_patterns[i]) could not be matched."))
        end

        numA = _to_int(inner_match[2])
        numB = _to_int(inner_match[4])

        hyperedge = Dict(A => numA, B => numB)
        push!(hyperedges, hyperedge)
    end

    return Label(hyperedges[1], hyperedges[2], hyperedges[3])
end

function Label(hyperedge::Hyperedge)
    return Label(copy(hyperedge), Dict(A => 0, B => 0), Dict(A => 0, B => 0))
end

function Label(hyperedge1::Hyperedge, hyperedge2::Hyperedge, int_state::State)
    @assert hyperedge1[int_state] > 0
    @assert hyperedge2[int_state] > 0

    left = copy(hyperedge1)
    left[int_state] -= 1
    right = copy(hyperedge2)
    right[int_state] -= 1
    int = Dict(A => 0, B => 0)
    int[int_state] = 1

    return Label(left, int, right)
end

function Base.show(io::IO, l::Label)
    if order(l) == 2
        return print(io,
                     "[ $(_to_str(l.left)) | $(_to_str(l.int)) | $(_to_str(l.right))]")
    else
        return print(io, "[ $(_to_str(l.left)) ]")
    end
end

function Base.getproperty(obj::Label, sym::Symbol)
    if sym === :left_total
        return Dict(A => obj.left[A] + obj.int[A], B => obj.left[B] + obj.int[B])
    elseif sym === :right_total
        return Dict(A => obj.right[A] + obj.int[A], B => obj.right[B] + obj.int[B])
    elseif sym === :int_state
        if obj.int[A] > 0 && obj.int[B] == 0
            return A
        elseif obj.int[B] > 0 && obj.int[A] == 0
            return B
        else
            return nothing
        end
    else
        return getfield(obj, sym)
    end
end

function order(l::Label)
    if l.int[A] != 0 || l.int[B] != 0 || l.right[A] != 0 || l.right[B] != 0
        return 2
    end
    if l.left[A] + l.left[B] > 1
        return 1
    end
    return 0
end

function Base.size(l::Label)
    if order(l) == 1
        return (l.left[A] + l.left[B],)
    end
    return (l.left[A] + l.left[B] + l.int[A] + l.int[B],
            l.right[A] + l.right[B] + l.int[A] + l.int[B])
end

let cache::Union{Vector{Label},Nothing} = nothing,
    cached_size::Union{Int64,Nothing} = nothing

    global function all_labels(max_size::Int64)
        if !isnothing(cache) && max_size == cached_size
            return cache
        end

        @assert max_size >= 2

        labels = Label[]

        # order zero
        push!(labels, Label("[A]"))
        push!(labels, Label("[B]"))

        # [An Bm | A / B | Ai Bj]
        for m in 0:max_size, n in 0:(max_size - m)
            if n + m < 2
                continue
            end
            # order one
            push!(labels, Label("[A$n B$m]"))

            # order two
            for j in 0:(max_size - 1), i in 0:(max_size - 1 - j)
                if i + j < 1
                    continue
                end
                # A in the intersection
                if n > 0
                    push!(labels, Label("[A$(n-1) B$m | A | A$i B$j]"))
                end
                # B in the intersection
                if m > 0
                    push!(labels, Label("[A$n B$(m-1) | B | A$i B$j]"))
                end
            end
        end

        labels = unique(labels)

        cache = labels
        cached_size = max_size

        return labels
    end
end

function _to_int(m)
    if m === nothing
        return 0
    elseif length(m) == 0
        return 1
    end
    return parse(Int64, m)
end

function _to_str(h::Hyperedge)
    A_str = ""
    B_str = ""

    if h[A] == 0
        A_str = ""
    elseif h[A] == 1
        A_str = "A"
    else
        A_str = "A$(h[A])"
    end

    if h[B] == 0
        B_str = ""
    elseif h[B] == 1
        B_str = "B"
    else
        B_str = "B$(h[B])"
    end

    return "$(A_str)$(B_str)"
end