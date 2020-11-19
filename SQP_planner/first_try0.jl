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



function sparse_jac_con(A_fx,B_fx,z)
    J = spzeros(nc,nz)
    for i = 1:T-1
        x = z[idx_x[i]]
        u = z[idx_u[i]]
        A = A_fx(x)
        B = B_fx(u)
        @show A
        @show B

        @infiltrate
        error()

        # J[idx_cons[i],idx_x[i]] = -A
        # J[idx_cons[i],idx_u[i]] = -B
        # J[idx_cons[i],idx_x[i+1]] = I(nx)
    end
    return J
end


function runit()
u = NaN*zeros(nu)
x = NaN*zeros(nx)

A_fx = x -> ForwardDiff.jacobian(x -> discrete_dynamics(x, u, dt),x)
B_fx = u -> ForwardDiff.jacobian(u -> discrete_dynamics(x, u, dt),u)

z = randn(nz)

@show z[idx_x[1]]
@show z[idx_u[1]]

J2 = sparse_jac_con(A_fx,B_fx,z)

return J2, J_c, z
end


J2,Jc,z = runit()
