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
struct ForwardEuler <: TimeIntegrator end

"""
    struct RK2

Second-order Runge-Kutta method.
"""
struct RK2 <: TimeIntegrator end
