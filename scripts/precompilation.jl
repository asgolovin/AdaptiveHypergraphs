"""
A file to precompile the most expensive packages like Makie. 
Used by PackageCompiler to create a sysimage with the packages already loaded.
"""

using GLMakie

fig = Figure()
display(fig)
ax = Axis(fig[1, 1])

# example from the mesh documentation 
vertices = [0.0 0.0;
            1.0 0.0;
            1.0 1.0;
            0.0 1.0]
faces = [1 2 3;
         3 4 1]
mesh!(ax, vertices, faces)

x = [1, 2, 3]
y = [1, 2, 3]
linesegments!(ax, x, y)
lines!(ax, x, y)
scatter!(ax, x, y)
hlines!(ax, 1.0)
empty!(fig)