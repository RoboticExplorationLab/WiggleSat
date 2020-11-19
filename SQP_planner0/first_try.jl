cd(dirname(@__DIR__))
Pkg.activate(".")
using OSQP, ForwardDiff, SuiteSparse, SparseArrays, LinearAlgebra
using Infiltrator

function discrete_dynamics(x,u,dt)
    px = x[1]
    py = x[2]
    θ  = x[3]

    dθ = u[1]
    v = u[2]
    #update the state
    return [px + v*dt*cos(θ);
            py + v*dt*sin(θ);
            θ  + dt*dθ]
end
nx = 3
nu = 2
T = 50
dt = 0.2
nz = (T)*nx + (T-1)*nu
nc = (T-1)*nx
idx_x = [(t-1)*(nx+nu) .+ (1:nx) for t = 1:T]
idx_u = [(t-1)*(nx+nu) .+ ((nx+1):(nx+nu)) for t = 1:T-1]
idx_cons = [(t-1)*(nx) .+ (1:nx) for t = 1:(T-1)]


function jac_test(A_fx,B_fx,z)

    for i = 1
        x = z[1:nx]
        u = z[(nx+1):(nx + nu)]
        @show A_fx(x)
        @show B_fx(u)
    end

end

function fwd_Test()
    u = NaN*zeros(nu)
    x = NaN*zeros(nx)

    A_fx = x -> ForwardDiff.jacobian(x -> discrete_dynamics(x, u, dt),x)
    B_fx = u -> ForwardDiff.jacobian(u -> discrete_dynamics(x, u, dt),u)

    # this way works?
    u = randn(nu)
    x = randn(nx)
    # this way just shows NaNs
    z = [randn(nu);randn(nx)]
    jac_test(A_fx,B_fx,copy(z))

    # this way works?
    u = randn(nu)
    x = randn(nx)
    jac_test(A_fx,B_fx,copy([x;u]))
end

fwd_Test()
