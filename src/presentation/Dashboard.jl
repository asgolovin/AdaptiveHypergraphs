abstract type Plot end

struct HGPlot <: Plot end
struct StateCountPlot <: Plot end

function display(hgplot::HGPlot)
end

function display(sc::StateCountPlot)
end

function create_figure()
    f = Figure(resolution = (600, 400))
    display(f)

end