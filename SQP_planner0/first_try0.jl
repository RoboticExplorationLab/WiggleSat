using ForwardDiff

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
    for i = 1

        # here i update what x and u are (real numbers, not Nan)
        x = z[idx_x[i]]
        u = z[idx_u[i]]

        # for some reason, these still show up as NaN
        @show A_fx(x)
        @show B_fx(u)



    end
end


function runit()

    # initially I declare these as NaN
    u = NaN*zeros(nu)
    x = NaN*zeros(nx)

    A_fx = x_2 -> ForwardDiff.jacobian(x_2 -> discrete_dynamics(x_2, u, dt),x_2)
    B_fx = u -> ForwardDiff.jacobian(u -> discrete_dynamics(x, u, dt),u)

    z = randn(nz)

    @show z[idx_x[1]]
    @show z[idx_u[1]]

    sparse_jac_con(A_fx,B_fx,z)

end

runit()
