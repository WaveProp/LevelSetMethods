"""
    abstract type TimeIntegrator
    
Methods for integrating and ordinary differential equation of the form `yₜ =
F(y,t)` forward in time.
"""
abstract type TimeIntegrator end

"""
    struct ForwardEuler <: TimeIntegrator
    
Forward Euler (explicit) integration method.
"""
Base.@kwdef struct ForwardEuler <: TimeIntegrator
end

"""
    struct RK2
    
Second-order Runge-Kutta method.         
"""
Base.@kwdef struct RK2 <: TimeIntegrator
end

"""
    struct RKLM2
    
Low-memory implementation of second-order Runge-Kutta method. In contrast to
`RK2`, the `RKLM2` method requires only one `buffer` for its time integration.
"""
Base.@kwdef struct RKLM2 <: TimeIntegrator
end



