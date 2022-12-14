export Hyperedge, Label

const Hyperedge = NamedTuple{(:A, :B)}

struct Label
    left::Hyperedge
    int::Hyperedge
    right::Hyperedge
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
    Label((A = 1, B = 1), (A = 0, B = 0), (A = 0, B = 0))

    julia> Label("[AB2 | B | A2B3]")
    Label((A = 1, B = 2), (A = 0, B = 1), (A = 2, B = 3))

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
                hyperedge = (A=0, B=0)
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

        hyperedge = (A=numA, B=numB)
        push!(hyperedges, hyperedge)
    end
    return Label(hyperedges[1], hyperedges[2], hyperedges[3])
end

function Base.show(io::IO, l::Label)
    if l.int.A != 0 || l.int.B != 0 || l.right.A != 0 || l.right.B != 0
        return print("""Label("[ $(_to_str(l.left)) | $(_to_str(l.int)) | $(_to_str(l.right)) ]")""")
    else
        return print("""Label("[ $(_to_str(l.left)) ]")""")
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

    if h.A == 0
        A_str = ""
    elseif h.A == 1
        A_str = "A"
    else
        A_str = "A$(h.A)"
    end

    if h.B == 0
        B_str = ""
    elseif h.B == 1
        B_str = "B"
    else
        B_str = "B$(h.B)"
    end

    return "$(A_str)$(B_str)"
end