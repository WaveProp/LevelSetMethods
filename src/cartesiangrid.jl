struct CartesianGrid{T}
    xrange::T
    yrange::T
end    

Base.size(g::CartesianGrid) = length(g.xrange), length(g.yrange)

gridsize(g::CartesianGrid) = step(g.xrange), step(g.yrange)
