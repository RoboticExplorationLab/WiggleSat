# cd(dirname(@__DIR__))
# Pkg.activate(".")
using OSQP, ForwardDiff, SuiteSparse, SparseArrays, LinearAlgebra
using Infiltrator
using BenchmarkTools
using StaticArrays
using Attitude, MATLAB
# function discrete_dynamics(x,u,dt)
#     px = x[1]
#     py = x[2]
#     θ  = x[3]
#
#     dθ = u[1]
#     v = u[2]
#     #update the state
#     return [px + v*dt*cos(θ);
#             py + v*dt*sin(θ);
#             θ  + dt*dθ]
# end
function discrete_dynamics(x,u,dt)
    px = x[1]
    py = x[2]
    vx = x[3]
    vy = x[4]

    ax = u[1]
    ay = u[2]
    #update the state
    return [vx;vy;ax;ay]
end





function constraint(z)

    c = zeros(eltype(z),(T-1)*nx)

    for i = 1:(T-1)
        xt   = z[idx_x[i]]
        ut   = z[idx_u[i]]
        xtp1 = z[idx_x[i+1]]

        c[idx_cons[i]] = xtp1 - discrete_dynamics(xt,ut,dt)
    end
    return c

end

function sparse_jac_con!(J,z)
    for i = 1:T-1
        # get current x and u
        x = SVector{nx}(z[idx_x[i]])
        u = SVector{nu}(z[idx_u[i]])

        # new closure every time, this seems bad
        A_fwd_fx(x) = discrete_dynamics(x,u,dt)
        B_fwd_fx(u) = discrete_dynamics(x,u,dt)

        # put them in the constraint jacobian
        J[idx_cons[i],idx_x[i]] = -ForwardDiff.jacobian(A_fwd_fx,x)
        J[idx_cons[i],idx_u[i]] = -ForwardDiff.jacobian(B_fwd_fx,u)
        J[idx_cons[i],idx_x[i+1]] = I(nx)
    end
    return J
end


nx = 4
nu = 2
T = 500
dt = 0.2
nz = (T)*nx + (T-1)*nu
nc = (T-1)*nx

idx_x = [(t-1)*(nx+nu) .+ (1:nx) for t = 1:T]
idx_u = [(t-1)*(nx+nu) .+ ((nx+1):(nx+nu)) for t = 1:T-1]
idx_cons = [(t-1)*(nx) .+ (1:nx) for t = 1:(T-1)]

function runit()


z = randn(nz)
J = spzeros(nc,nz)
# @time ForwardDiff.jacobian(constraint,z)
sparse_jac_con!(J,z)
# @infiltrate
# @show norm(Jc-J)

# initial conditions
x0 = [20;30;0;0]
# create OSQP problem
z = [x0;0.001*randn(nz-nx)]

Q = sparse(10*I(nx))
Qf = copy(Q)
R = sparse(10*I(nu))
P = blockdiag(kron(sparse(I(T-1)),blockdiag(Q,R)),Qf)
q = zeros(nz)
A = spzeros(nc,nz)

# get jacobian
sparse_jac_con!(A,z)

# equality
A_lower = spzeros(nx,nz)
A_lower[1:nx,1:nx] = sparse(I(nx))
upper = (A*z - constraint(z))
lower = copy(upper)
A = [A; A_lower]
upper = [upper;x0]
lower = [lower;x0]
# osqp stuff
m = OSQP.Model()

OSQP.setup!(m; P = P, q=q, A = A, l = lower, u = upper,eps_abs = 1e-8,eps_rel = 1e-8)

for i = 1:1

    results = OSQP.solve!(m)
    z=copy(results.x)

    # generate new stuff
    A = spzeros(nc,nz)
    # get jacobian
    sparse_jac_con!(A,z)
    # equality
    A_lower = spzeros(nx,nz)
    A_lower[1:nx,1:nx] = sparse(I(nx))
    upper = (A*z - constraint(z))
    lower = copy(upper)
    A = [A; A_lower]
    upper = [upper;x0]
    lower = [lower;x0]

    # @infiltrate
    # error()
    OSQP.update!(m,l = lower, u = upper)
    OSQP.update!(m, Ax = A.nzval, Ax_idx = (collect( 1: length(A.nzval))))

end


X = fill(zeros(nx),T)
U = fill(zeros(nu),T-1)
for i = 1:T-1
    X[i] = z[idx_x[i]]
    U[i] = z[idx_u[i]]
end
X[T] = z[idx_x[T]]


Xm = mat_from_vec(X)
Um = mat_from_vec(U)

mat"
figure
hold on
plot($Xm')
hold off
"
mat"
figure
hold on
plot($Um')
hold off
"
mat"
figure
hold on
plot($Xm(1,:),$Xm(2,:))
hold off
"


@show norm(constraint(z))


end


runit()


# function setup_update_matrices()
#     options = Dict(:verbose => false,
#                    :eps_abs => 1e-08,
#                    :eps_rel => 1e-08,
#                    :polish => false,
#                    :check_termination => 1)
#
#     # seed!(1)
#
#     n = 5
#     m = 8
#     p = 0.7
#     Pt = sprandn(n, n, p)
#     Pt_new = copy(Pt)
#     P = Pt * Pt' + sparse(I, n, n)
#     #  PtI = findall(!iszero, P)
#     #  (Pti, Ptj) = (getindex.(PtI, 1), getindex.(PtI, 2))
#     Ptx = copy(Pt.nzval)
#
#     Pt_newx = Ptx + 0.1 * randn(length(Ptx))
#     # Pt_new = sparse(Pi, Pj, Pt_newx)
#     P_new = Pt_new * Pt_new' + sparse(I, n, n)
#     q = randn(n)
#     A = sprandn(m, n, p)
#
#     (Ai, Aj, _) = findnz(A)
#     Ax = copy(A.nzval)
#     A_newx = Ax + randn(length(Ax))
#     A_new = sparse(Ai, Aj, A_newx)
#
#     # A_new = copy(A)
#     # A_new.nzval += randn(length(A_new.nzval))
#     l = zeros(m)
#     u = 30 .+ randn(m)
#
#
#     problem = Dict()
#     problem[:P] = P
#     problem[:P_new] = P_new
#     problem[:q] = q
#     problem[:A] = A
#     problem[:A_new] = A_new
#     problem[:l] = l
#     problem[:u] = u
#     problem[:m] = m
#     problem[:n] = n
#     return problem, options
# end
#
# problem, options = setup_update_matrices()
#
# (n, m, P, q, A, l, u) = (problem[:n], problem[:m],
#                              problem[:P], problem[:q], problem[:A], problem[:l], problem[:u])
# A_new = problem[:A_new]
#
# model = OSQP.Model()
# OSQP.setup!(model; P=P, q=q, A=A, l=l, u=u, options...)
#
#                 # Update matrix
# A_new_idx = collect(1 : length(A_new.nzval))
# X = fill(zeros(nx),T)
# U = fill(zeros(nu),T-1)
# X[1] = randn(nx)
#
# for i = 1:T-1
#     U[i] = randn(nu)
#     X[i+1] = discrete_dynamics(X[i],U[i],dt)
# end
#
# z = zeros(nz)
# for i = 1:T-1
#     z[idx_x[i]] = X[i]
#     z[idx_u[i]] = U[i]
# end
# z[idx_x[T]] = X[T]

















# test constraint function
# X = fill(zeros(nx),T)
# U = fill(zeros(nu),T-1)
# X[1] = randn(nx)
#
# for i = 1:T-1
#     U[i] = randn(nu)
#     X[i+1] = discrete_dynamics(X[i],U[i],dt)
# end
#
# z = zeros(nz)
# for i = 1:T-1
#     z[idx_x[i]] = X[i]
#     z[idx_u[i]] = U[i]
# end
# z[idx_x[T]] = X[T]
