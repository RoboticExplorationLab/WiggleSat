using JuMP
using LinearAlgebra
# using OSQP
using Mosek
using Attitude
using MosekTools
# using COSMO
using JLD2
@load "/Users/kevintracy/devel/WiggleSat/convex_planner/orbit_data_2.jld2" τ_hist B_hist_b J

N = length(τ_hist)
dt = 10.0
# tol = 1e-9
model = Model(Mosek.Optimizer)
set_optimizer_attribute(model, "MSK_DPAR_INTPNT_CO_TOL_DFEAS",1e-15)
set_optimizer_attribute(model, "MSK_DPAR_INTPNT_CO_TOL_REL_GAP",1e-15)


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


t_vec = (0:dt:(N-1)*dt)/(3600)
t_vec_slim = (0:dt:(N-2)*dt)/3600
using MATLAB
mat"
figure
hold on
title('Arm Angles')
plot($t_vec,rad2deg($θ'))
ylabel('Arm Angle (deg)')
xlabel('Time (hours)')
hold off
"
mat"
figure
hold on
title('Arm Angle Derivative')
plot($θ̇')
hold off
"
# mat"
# figure
# hold on
# title('Arm Torque','FontSize', 18)
# plot($t_vec_slim,$J_arms*$α','linewidth',3)
# h=fill([1,1,1.57,1.57],[-0.0025,0.004,0.004,-0.0025],'b');
#     h.FaceAlpha=0.2;
#     h.EdgeAlpha = 0.0;
#     h=fill([2.58,2.58,3.1667,3.1667],[-0.0025,0.004,0.004,-0.0025],'b');
#         h.FaceAlpha=0.2;
#         h.EdgeAlpha = 0.0;
# xlabel('Time (hours)')
# xlim([0 $t_vec_slim(end)])
# legend([h(1)],'Eclipse','FontSize',15)
# ylim([min(vec($J_arms*$α)) max(vec($J_arms*$α))])
# ylabel('Arm Torque (Nm)')
# pbaspect([2 1 1])
# hold off
# saveas(gcf,'armtorque.eps','epsc')
# "
# mat"
# figure
# hold on
# title('Magnetic Moment','FontSize', 18)
# plot($t_vec_slim,$m','linewidth',3)
# h=fill([1,1,1.57,1.57],[-0.0025,0.004,0.004,-0.0025],'b');
#     h.FaceAlpha=0.2;
#     h.EdgeAlpha = 0.0;
#     h=fill([2.58,2.58,3.1667,3.1667],[-0.0025,0.004,0.004,-0.0025],'b');
#         h.FaceAlpha=0.2;
#         h.EdgeAlpha = 0.0;
# xlabel('Time (hours)')
# xlim([0 $t_vec(end)])
# legend([h(1)],'Eclipse','FontSize',15)
# ylim([min(vec($m)) max(vec($m))])
# ylabel('Magnetic Moment (A/m^2)')
# pbaspect([2 1 1])
# hold off
# saveas(gcf,'magnet.eps','epsc')
# "
mat"
figure
hold on
subplot(2,1,1)
hold on
title('Arm Torque','FontSize', 15)
plot($t_vec_slim,$J_arms*$α','linewidth',3)
h=fill([1,1,1.57,1.57],[-0.0025,0.004,0.004,-0.0025],'b');
    h.FaceAlpha=0.2;
    h.EdgeAlpha = 0.0;
    h=fill([2.58,2.58,3.1667,3.1667],[-0.0025,0.004,0.004,-0.0025],'b');
        h.FaceAlpha=0.2;
        h.EdgeAlpha = 0.0;
xlabel('Time (hours)')
xlim([0 $t_vec_slim(end)])
legend([h(1)],'Eclipse','FontSize',15)
ylim([min(vec($J_arms*$α)) max(vec($J_arms*$α))])
ylabel('Arm Torque (Nm)')
%pbaspect([2 1 1])
hold off
%saveas(gcf,'armtorque.eps','epsc')

subplot(2,1,2)
hold on
title('Magnetic Moment','FontSize', 15)
plot($t_vec_slim,$m','linewidth',3)
h=fill([1,1,1.57,1.57],[-0.0025,0.004,0.004,-0.0025],'b');
    h.FaceAlpha=0.2;
    h.EdgeAlpha = 0.0;
    h=fill([2.58,2.58,3.1667,3.1667],[-0.0025,0.004,0.004,-0.0025],'b');
        h.FaceAlpha=0.2;
        h.EdgeAlpha = 0.0;
xlabel('Time (hours)')
xlim([0 $t_vec(end)])
legend([h(1)],'Eclipse','FontSize',15)
ylim([min(vec($m)) max(vec($m))])
ylabel('Magnetic Moment (A/m^2)')
%pbaspect([2 1 1])
hold off
saveas(gcf,'convex_slew_data.eps','epsc')
"

tau = mat_from_vec(τ_hist)

mat"
figure
hold on
title('Disturbance Torques','FontSize', 18)
plot($t_vec,$tau','linewidth',3)
ylabel('Disturbance Torques (Nm)')
xlabel('Time (hours)')
xlim([0 $t_vec(end)])
hold off
saveas(gcf,'torques.png')
"

println("test")
# torque_matching_errors = zeros(N-1)
# for i = 1:N-1
# torque_matching_errors[i] = norm((skew_from_vec(B_hist_b[i])*m[:,i] + J_arms*α[:,i]) -τ_hist[i])
# end
#
# mat"
# figure
# hold on
# plot($torque_matching_errors)
# hold off
# "
#
# θm = vec_from_mat(θ)

# @save "state_hist_for_vis.jld2" θm

## now let's try a real sim
# using Dierckx
#
#
# # spl = Spline1D(x, y)
#
# τ_d1 = Spline1D(0:dt:dt*(N-1),mat_from_vec(τ_hist)[1,:])
# τ_d2 = Spline1D(0:dt:dt*(N-1),mat_from_vec(τ_hist)[2,:])
# τ_d3 = Spline1D(0:dt:dt*(N-1),mat_from_vec(τ_hist)[3,:])
# B_eci1 = Spline1D(0:dt:dt*(N-1),mat_from_vec(B_hist_b)[1,:])
# B_eci2 = Spline1D(0:dt:dt*(N-1),mat_from_vec(B_hist_b)[2,:])
# B_eci3 = Spline1D(0:dt:dt*(N-1),mat_from_vec(B_hist_b)[3,:])
#
# function B_eci(t)
#     return [B_eci1(t);B_eci2(t);B_eci3(t)]
# end
# function τ_dist(t)
#     return [τ_d1(t);τ_d2(t);τ_d3(t)]
# end
#
# nx = 6
# nu = 6
#
# dt = 10.0
#
# X_plan = [[zeros(6);θ[:,i];θ̇[:,i]] for i = 1:N]
# U_plan = [[m[:,i];α[:,i]] for i = 1:(N-1)]
# function rk4(f, u, x_n, h,t_n)
#     """Vanilla RK4"""
#     k1 = h*f(x_n,u,t_n)
#     k2 = h*f(x_n+k1/2, u,t_n + h/2)
#     k3 = h*f(x_n+k2/2, u, t_n + h/2)
#     k4 = h*f(x_n+k3, u, t_n + h)
#     return (x_n + (1/6)*(k1+2*k2+2*k3 + k4))
# end
# function dynamics(x,u,t)
#     ᴺpᴮ = x[1:3]
#     ω = x[4:6]
#
#     m = u[1:3]
#     α = u[4:6]
#
#     B_b = dcm_from_p(ᴺpᴮ)'*B_eci(t)
#     τ_d = τ_dist(t)
#     τ = τ_d - skew_from_vec(B_b)*m - J_arms*α
#     # τ = -J_arms*α
#     # @infiltrate
#     # error()
#     return [pdot_from_w(ᴺpᴮ,ω);invJ*(τ - ω × (J*ω))]
# end
#
# X = fill(zeros(nx),N)
# # X[1][1:3] = 0.001*randn(3)
# invJ = inv(J)
# nn = 600
# for i = 1:nn
#     t = (i-1)*dt
#
#     δx = X[i][1:6] - X_plan[i][1:6]
#     kp = 1e-4
#     kd = 20e-4
#     u_fb = kp*δx[1:3] + kd*δx[4:6]
#
#     u = U_plan[i] + [zeros(3);u_fb]
#     # u = [zeros(3);u_fb]*1e-4
#
#     X[i+1] = rk4(dynamics,u,X[i],dt,t)
# end
#
# Xm = mat_from_vec(X)
#
# mat"
# figure
# hold on
# plot($Xm(1:3,1:$nn)')
# hold off
# "
