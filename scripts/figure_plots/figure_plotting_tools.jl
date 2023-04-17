using CairoMakie
using CSV
using DataFrames
using Colors
using ColorSchemes
using Polynomials
using DrWatson
using REPL.TerminalMenus

TUM_primary_blue = colorant"#0065BD"
TUM_secondary_light_blue = colorant"#005293"
TUM_secondary_dark_blue = colorant"#003359"
TUM_secondary_dark_gray = colorant"#333333"
TUM_secondary_medium_gray = colorant"#808080"
TUM_secondary_light_gray = colorant"#CCCCC6"
TUM_accent_gray = colorant"#DAD7CB"
TUM_accent_orange = colorant"#E37222"
TUM_accent_green = colorant"#A2AD00"
TUM_accent_very_light_blue = colorant"#98C6EA"
TUM_accent_light_blue = colorant"#64A0C8"

tum_scheme = ColorScheme([TUM_accent_orange,
                          TUM_accent_green,
                          TUM_accent_light_blue,
                          TUM_secondary_light_gray])

screen_config = Dict(:pt_per_unit => 1)
CairoMakie.activate!(; screen_config...)

mytheme = Theme(;
                Axis=(titlesize=12,
                      xticksize=3,
                      yticksize=3,
                      ylabelpadding=2,
                      xticklabelsize=10,
                      yticklabelsize=10,
                      xticklabelspace=10.0,
                      xticklabelsvisible=true,
                      yticklabelsvisible=true,
                      xlabelsize=12,
                      ylabelsize=12))
set_theme!(mytheme)

function fit_parabola(x, y, fixed_at_one::Bool=false)
    if fixed_at_one
        # p(x) = a x^2 - a
        X = x .^ 2 .- 1
        coeff = (X' * X)^(-1) * X' * y
        sigmasq = 1 / (length(y) - 2) * (y - X * coeff)' * (y - X * coeff)
        cov = (X' * X)^(-1) * sigmasq
        p = Polynomial([-coeff, 0.0, coeff])
    else
        # p(x) = c x^2 + a
        X = ones((length(x), 2))
        X[1:end, 2] = x .^ 2
        coeff = (X' * X)^(-1) * X' * y
        sigmasq = 1 / (length(y) - 2) * (y - X * coeff)' * (y - X * coeff)
        cov = (X' * X)^(-1) * sigmasq
        p = Polynomial([coeff[1], 0.0, coeff[2]])
    end
    return p, cov
end

struct DataFolder
    folder::String
    type::Symbol # :simulation or :batch
end

function Base.iterate(f::DataFolder)
    folders = []
    for folder in readdir(f.folder; join=true)
        if !isdir(folder)
            continue
        end
        if f.type == :simulation
            push!(folders, DataFolder(folder, :batch))
        else
            push!(folders, folder)
        end
    end

    if length(folders) == 0
        return nothing
    end

    if f.type == :simulation
        folder_names = [splitpath(f.folder)[end] for f in folders]
        ids = [parse(Int64, match(r"batch_([0-9]+)", name)[1])
               for name in folder_names]
    elseif f.type == :batch
        folder_names = [splitpath(f)[end] for f in folders]
        ids = [parse(Int64, match(r"run_([0-9]+)", name)[1])
               for name in folder_names]
    end

    result = (folders[1], ids[1])
    state = (folders[2:end], ids[2:end])

    return (result, state)
end

function Base.iterate(f::DataFolder, state)
    folders, ids = state

    if length(folders) == 0
        return nothing
    end

    result = (folders[1], ids[1])
    state = (folders[2:end], ids[2:end])

    return (result, state)
end

function prompt_for_save(filename, fig; pt_per_unit=0.666)
    ans = Base.prompt("Press space and then twice Enter")
    options = ["no", "yes"]
    menu = RadioMenu(options)

    choice = request("Save the figure to file $filename?", menu)
    if choice == 2
        save(filename, fig; pt_per_unit)
    end
end

function load_data(filename, run_folder)
    path = joinpath(run_folder, filename)
    df = CSV.read(path, DataFrame; delim=", ")
    time = df.index
    value = df.value
    return (time, value)
end

"""
Compute the maximum of two vectors element-wise, works for 
vectors of different length. 

# Example: 
    julia> safe_max([1, 2, 3, 4], [43, 77, -1])
    [43, 77, 3, 4]
"""
function safe_max(vector1, vector2)
    l1 = length(vector1)
    l2 = length(vector2)
    if l1 == l2
        return max.(vector1, vector2)
    elseif l1 < l2
        return vcat(max.(vector1, vector2[1:l1]), vector2[(l1 + 1):end])
    else
        return vcat(max.(vector1[1:l2], vector2), vector1[(l2 + 1):end])
    end
end

"""
Compute the minimum of two vectors element-wise, works for 
vectors of different length. 

# Example: 
    julia> safe_min([1, 2, 3, 4], [43, 77, -1])
    [1, 2, -1, 4]
"""
function safe_min(vector1, vector2)
    l1 = length(vector1)
    l2 = length(vector2)
    if l1 == l2
        return min.(vector1, vector2)
    elseif l1 < l2
        return vcat(min.(vector1, vector2[1:l1]), vector2[(l1 + 1):end])
    else
        return vcat(min.(vector1[1:l2], vector2), vector1[(l2 + 1):end])
    end
end