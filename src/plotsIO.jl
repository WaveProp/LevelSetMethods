# recipes for Plots
@recipe function f(eq::LevelSetEquation)
    ϕ = levelset(eq)
    t = current_time(eq)
    str = @sprintf "t = %.2f" t
    @series begin
        title --> str
        ϕ
    end
end


@recipe function f(ls::LevelSet)
    ϕ = levelset_function(ls)
    s = levelset_sign(ls)
    N = ambient_dimension(ls)
    if N == 2 # 2d contour plot
        if s == 0
            seriestype --> :contour
        else
            seriestype --> :contourf
        end
        levels --> [0,0]
        aspect_ratio --> :equal
        colorbar --> false
        # seriescolor --> :black
        m = vals_mesh(ϕ)
        # Note: the vals of ϕ need be transposed because contour expects the
        # matrix to have rows representing the x vals and columns expecting
        # the y value.
        xgrid = grids(m,1)
        ygrid = grids(m,2)
        return xgrid,ygrid,-transpose(vals(ϕ))
    else
        notimplemented()
    end
end
