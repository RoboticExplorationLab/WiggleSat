using LinearAlgebra, ForwardDiff, COSMO, Convex



m = 4
n = 20

A = randn(m,n)
b = abs.(randn(m))

x_cvx = Variable(n)

Q = randn(n,n);Q = Q'*Q;
q = randn(n)
cons = Constraint[ A*x_cvx - b <= 0]

problem = minimize(.5*quadform(x_cvx,Q) + dot(q,x_cvx),cons)


solve!(problem, () -> COSMO.Optimizer(eps_abs = 1e-8,eps_rel = 1e-8))


# now let's try with barrier

function c(x)
    return A*x - b
end
function f(x)
    return 0.5*x'*Q*x + dot(q,x)
end
function L(x,λ,I_mu)
    c1 = c(x)
    return f(x) + dot(λ,c1) + .5*c1'*I_mu*c1
end
function get_I(x,λ,μ)
    I_mu = float(I(m))
    c1 = c(x)
    for i = 1:m
        if (c1[i] < 0) && (λ[i] == 0)
            I_mu[i,i] = 0
        else
            I_mu[i,i] = μ[i]
        end
    end
    return I_mu
end
function λ_update!(λ,x,μ)
    c1 = c(x)
    for i = 1:m
        # @infiltrate
        # error()
        λ[i] = max(0,λ[i] + μ[i]*c1[i])
    end
end
function aug()

    x = zeros(n)
    λ = zeros(m)
    μ = ones(m)

    ϕ = 5
    # μ = 1

    for i = 1:100

        I_mu = get_I(x,λ,μ)

        _L(x) = L(x,λ,I_mu)

        H = ForwardDiff.hessian(_L,x)
        g = ForwardDiff.gradient(_L,x)

        #-------------------------------------------
        # this is where the problem is solved (OSQP)
        v = - H\g
        x = x +v
        #-------------------------------------------
        λ_update!(λ,x,μ)

        μ *= ϕ

        @show norm(v)
        if norm(v) < 1e-6
            break
        end


    end
    return x
end

x = aug()
