using Convex, MosekTools, Mosek, Attitude, COSMO


@load "/Users/kevintracy/devel/WiggleSat/convex_planner/orbit_data.jld2" τ_hist B_hist_b J


G = 1e7

J_arm = 0.125

max_moments = .01*ones(3)

N = length(τ_hist)

# states
ϕ = Variable(3,N)
ω = Variable(3,N)
θ = Variable(3,N)
θ_dot = Variable(3,N)

# controls
m = Variable(3,N-1)
α = Variable(3,N-1)

dt = 10.0

cons = Constraint[ θ[:,1] == zeros(3)]
push!(cons, ϕ[:,1] == zeros(3))
push!(cons, ω[:,1] == zeros(3))
push!(cons, θ_dot[:,1] == zeros(3))

for i = 1:N-1
    # arm dynamics
    push!(cons, θ_dot[:,i+1] == θ_dot[:,i] + dt*α[:,i])

    # arm kinematics
    push!(cons, θ[:,i+1] == θ[:,i] + dt*θ_dot[:,i])

    # spacecraft dynamics
    push!(cons, ω[:,i+1] == ω[:,i] + dt*G*(τ_hist[i] -
                skew_from_vec(B_hist_b[i])*m[:,i] + J_arm*α[:,i] ))

    # spacecraft kinematics
    push!(cons, ϕ[:,i+1] == ϕ[:,i] + dt*ω[:,i])
end

for i = 356:569
    push!(cons, m[:,i] == zeros(3))
end

# add magnetic moments constraints
for i = 1:N-1
    push!(cons, m <= 0.05)
    push!(cons, m >= -0.05)
end

prob = minimize( sumsquares(ϕ) + sumsquares(ω) + sumsquares(θ) +
                 sumsquares(θ_dot) + sumsquares(m) + sumsquares(α) ,cons)


solve!(prob, () -> Mosek.Optimizer())
# solve!(prob, () -> COSMO.Optimizer(max_iter = 20000))


# states
ϕ = evaluate(ϕ)
ω = evaluate(ω)
θ = evaluate(θ)
θ_dot = evaluate(θ_dot)

# controls
m = evaluate(m)
α = evaluate(α)

mat"
figure
hold on
title('Attitude')
plot($ϕ')
hold off
"
mat"
figure
hold on
title('Angular Velocity')
plot($ω')
hold off
"
mat"
figure
hold on
title('Arm Angles')
plot($θ')
hold off
"
mat"
figure
hold on
title('Arm Angle Derivative')
plot($θ_dot')
hold off
"
mat"
figure
hold on
title('Arm Angle Acceleration')
plot($α')
hold off
"
mat"
figure
hold on
title('Magnetic Moment')
plot($m')
hold off
"
