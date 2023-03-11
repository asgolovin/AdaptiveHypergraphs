export AbstractMotif, OrderZeroMotif, OrderOneMotif, OrderTwoMotif, State, order, size,
       all_motifs

@enum State::Bool A = false B = true

const MotifPart = NamedTuple{(:A, :B),Tuple{Int64,Int64}}

abstract type AbstractMotif end

struct OrderZeroMotif <: AbstractMotif
    state::State
end

struct OrderOneMotif <: AbstractMotif
    A::Int64
    B::Int64

    function OrderOneMotif(a, b)
        @assert a >= 0 && b >= 0
        @assert a + b >= 2
        return new(a, b)
    end
end

function Base.size(motif::OrderOneMotif)
    return motif.A + motif.B
end

"""
    ind(motif::OrderZeroMotif)

Return the index of the motif. 

The motifs of order zero and one are enumerated like [A], [B], [A2], [AB], [B2], [A3], [A2B], [AB2], [B3], ...
"""
function index(motif::OrderZeroMotif)
    return motif.state == A ? 1 : 2
end

"""
    index(motif::OrderOneMotif)

Return the index of the motif. 

The motifs of order zero and one are ordered like [A], [B], [A2], [AB], [B2], [A3], [A2B], [AB2], [B3], ...
"""
function index(motif::OrderOneMotif)
    motif_size = size(motif)
    # the index at which the motifs of size `motif_size` start
    offset = sum([s for s in 2:motif_size]) + 1
    return offset + motif.B
end

function is_active(motif::OrderOneMotif)
    return motif.A > 0 && motif.B > 0
end

struct OrderTwoMotif <: AbstractMotif
    left::MotifPart
    int::MotifPart
    right::MotifPart

    function OrderTwoMotif(left::MotifPart, int::MotifPart, right::MotifPart)
        @assert left.A + int.A + left.B + int.B >= 2
        @assert int.A + int.B >= 1
        @assert int.A + right.A + int.B + right.B >= 1
        @assert left.A >= 0 && left.B >= 0
        @assert int.A >= 0 && int.B >= 0
        @assert right.A >= 0 && right.B >= 0

        # swap the left and right parts around such that the hyperedge with the 
        # larger index comes first
        left_motif = OrderOneMotif(left.A + int.A, left.B + int.B)
        right_motif = OrderOneMotif(int.A + right.A, int.B + right.B)

        if index(right_motif) > index(left_motif)
            right, left = left, right
        end

        return new(left, int, right)
    end
end

function OrderTwoMotif(left::Tuple{Int64,Int64},
                       int::Tuple{Int64,Int64},
                       right::Tuple{Int64,Int64})
    return OrderTwoMotif(MotifPart(left), MotifPart(int), MotifPart(right))
end

function Base.show(io::IO, m::OrderZeroMotif)
    return print(io, "[ $(m.state) ]")
end

function Base.show(io::IO, m::OrderOneMotif)
    motif_part = MotifPart((A=m.A, B=m.B))
    return print(io, "[ $(_to_str(motif_part)) ]")
end

function Base.show(io::IO, m::OrderTwoMotif)
    left_str = _to_str(m.left)
    int_str = _to_str(m.int)
    right_str = _to_str(m.right)
    return print(io, "[ $left_str | $int_str | $right_str ]")
end

function Base.getproperty(obj::OrderTwoMotif, sym::Symbol)
    if sym === :left_motif
        return OrderOneMotif(obj.left.A + obj.int.A, obj.left.B + obj.int.B)
    elseif sym === :right_motif
        return OrderOneMotif(obj.right.A + obj.int.A, obj.right.B + obj.int.B)
    else
        return getfield(obj, sym)
    end
end

order(m::OrderZeroMotif) = 0
order(m::OrderOneMotif) = 1
order(m::OrderTwoMotif) = 2

let cache::Union{Vector{AbstractMotif},Nothing} = nothing,
    cached_size::Union{Int64,Nothing} = nothing,
    cached_int_size::Union{Int64,Nothing} = nothing

    global function all_motifs(max_size::Int64, int_size::Int64=max_size - 1)
        if !isnothing(cache) && max_size == cached_size && int_size == cached_int_size
            return cache
        end

        @assert max_size >= 2
        @assert 1 <= int_size < max_size

        motifs = AbstractMotif[]

        # order zero
        push!(motifs, OrderZeroMotif(A))
        push!(motifs, OrderZeroMotif(B))

        # order one
        for size in 2:max_size, b in 0:size
            a = size - b
            push!(motifs, OrderOneMotif(a, b))
        end

        # order two
        # [An Bm | Ar Bs | Aa Bb]
        # n + m + r + s = size_left, r + s + a + b = size_right
        # m + s = B_left, s + b = B_right
        # r ∈ [0, min(A_left, A_right)], s ∈ [0, min(B_left, B_right)]
        # 1 <= r + s <= int_size

        # Iterate over the sizes of the two hyperedges first
        for size_left in 2:max_size, size_right in 2:size_left
            # Then over the composition of A- and B-nodes in every hyperedge
            for B_left in 0:size_left, B_right in 0:size_right
                A_left = size_left - B_left
                A_right = size_right - B_right
                # Then over the possible ways how those hyperedges can intersect
                for r in 0:min(A_left, A_right), s in 0:min(B_left, B_right)
                    # The intersection has to be non-empty and not exceed int_size
                    if !(1 <= r + s <= int_size)
                        continue
                    end
                    n = A_left - r
                    m = B_left - s
                    a = A_right - r
                    b = B_right - s
                    push!(motifs, OrderTwoMotif((n, m), (r, s), (a, b)))
                end
            end
        end

        motifs = unique(motifs)

        cache = motifs
        cached_size = max_size
        cached_int_size = int_size

        return motifs
    end
end

function _to_str(motif_part::MotifPart)
    A_str = ""
    B_str = ""

    if motif_part.A == 0
        A_str = ""
    elseif motif_part.A == 1
        A_str = "A"
    else
        A_str = "A$(motif_part.A)"
    end

    if motif_part.B == 0
        B_str = ""
    elseif motif_part.B == 1
        B_str = "B"
    else
        B_str = "B$(motif_part.B)"
    end

    return "$(A_str)$(B_str)"
end

function motif_from_string(motif_string)
    # TODO
end