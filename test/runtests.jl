using AdaptiveHypergraphs
using Test

for file in filter(f -> endswith(f, ".jl"), readdir(@__DIR__))
    if file == "runtests.jl"
        continue
    end

    test_name = replace(file, "_" => " ", ".jl" => "")
    @testset "$test_name" begin
        t = time()
        include(file)
        println("$(test_name) took $(round(time() - t; digits = 1)) seconds.")
    end
end