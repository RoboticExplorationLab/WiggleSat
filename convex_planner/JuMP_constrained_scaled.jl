using JuMP
using OSQP
using Mosek
using MosekTools
using COSMO
@load "/Users/kevintracy/devel/WiggleSat/convex_planner/orbit_data_2.jld2" τ_hist B_hist_b J

N = length(τ_hist)
dt = 10.0
tol = 1e-9
model = Model(Mosek.Optimizer)
set_optimizer_attribute(model, "MSK_DPAR_INTPNT_CO_TOL_DFEAS",1e-12)
set_optimizer_attribute(model, "MSK_DPAR_INTPNT_CO_TOL_REL_GAP",1e-12)

# model = Model(COSMO.Optimizer)
# set_optimizer_attribute(model, "max_iter",20000)
# model = Model(OSQP.Optimizer)
# set_optimizer_attribute(model, "eps_abs", tol)
# set_optimizer_attribute(model, "eps_rel", tol)
# set_optimizer_attribute(model, "eps_prim_inf", tol)
# set_optimizer_attribute(model, "eps_dual_inf", tol)
# set_optimizer_attribute(model, "max_iter", 20000)
# set_optimizer_attribute(model, "polish", 1)

# create variables
@variable(model, α[1:3,1:N-1])
@variable(model, m[1:3,1:N-1])
@variable(model, θ[1:3,1:N])
@variable(model, θ̇[1:3,1:N])

# create torque matching constraint
J_arms = 0.03
sf = 1e6
for i = 1:N-1
    @constraint(model, sf*(skew_from_vec(B_hist_b[i])*m[:,i] + J_arms*α[:,i]) .== sf*τ_hist[i])
end

# arm stops
@constraint(model, θ .<=  float(pi))
@constraint(model, θ .>= -float(pi))

# make sure there is no mag torque during eclipse
for i = 1:N-1
    if (i >= 356 && i<=569) || (i >= 930 && i<=1141)
        @constraint(model, m[:,i] .== 0)
    else
        @constraint(model, m[:,i] .<=  0.01)
        @constraint(model, m[:,i] .>= -0.01)
    end
end

# dynamics of the arm
A = Array([zeros(3,3) I(3);zeros(3,6)])
B = Array([zeros(3,3);I(3)])
exp_discrete = exp([A B; zeros(3,9)]*dt)
Ad = exp_discrete[1:6,1:6]
Bd = exp_discrete[1:6,7:9]
for i = 1:N-1
    @constraint(model, [θ[:,i+1];θ̇[:,i+1]] .== Ad*[θ[:,i];θ̇[:,i]] + Bd*α[:,i] )
end

# initial conditions
@constraint(model, θ[:,1] .== 0)
@constraint(model, θ̇[:,1] .== 0)

# objective expression
θ_weight = 1e-10
θ̇_weight = 1e-10
α_weight = 1e3
m_weight = 1e-7
objective_exp = @expression(model, θ_weight*(θ[:,N])'*(θ[:,N]) +
                                  θ̇_weight *(θ̇[:,N])'*(θ̇[:,N]))
for i = 1:N-1
    add_to_expression!(objective_exp, θ_weight*(θ[:,i])'*(θ[:,i]) +
                                      θ̇_weight*(θ̇[:,i])'*(θ̇[:,i]) +
                                      α_weight*(α[:,i])'*(α[:,i]) +
                                      m_weight*(m[:,i])'*(m[:,i]) )
end
@objective(model, Min, objective_exp)

optimize!(model)


θ = value.(θ)
θ̇ = value.(θ̇)
m = value.(m)
α = value.(α)


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


torque_matching_errors = zeros(N-1)
for i = 1:N-1
torque_matching_errors[i] = norm((skew_from_vec(B_hist_b[i])*m[:,i]/sf + J_arms*α[:,i]/sf) -τ_hist[i])
end

mat"
figure
hold on
plot($torque_matching_errors)
hold off
"
