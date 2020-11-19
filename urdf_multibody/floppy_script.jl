using Pkg, LinearAlgebra

cd("/Users/kevintracy/julia_research/urdf_multibody/RBD_env")
Pkg.activate(".")

using RigidBodyDynamics, RigidBodySim
using MeshCat, MeshCatMechanisms

urdf = "/Users/kevintracy/julia_research/urdf_multibody/floppy_sc.xml"

# function I_from_abc(a,b,c,m)
#
#     ixx = (1/12)*m*(b^2+c^2)
#     iyy = (1/12)*m*(a^2+c^2)
#     izz = (1/12)*m*(a^2+b^2)
#
#     return diagm([ixx;iyy;izz])
# end
#
# I_bus = I_from_abc(3,1,1,10.0)
#
# I_sc = I_from_abc(2,3,.1,1.0)
# rob = parse_urdf(urdf,floating = true)
rob = parse_urdf(urdf,gravity = [0.0,0.0,0.0],floating = true)

state = MechanismState(rob)

# theta0 = deg2rad(10)
# r0 = normalize(randn(3))
# q0 = [cos(theta0/2);r0*sin(theta0/2)]
set_configuration!(state,[1;0;0;0;0;0;0;.1;-.5;.3])
set_velocity!(state,[0;0;0;0;0;0;0;0;0])
# set_velocity!(state,[0;0;0;0;.0;.0;0.0])
# set_velocity!(state,[.05])


# result = DynamicsResult(rob)
# #
# xdot = zeros(19)
# a = dynamics!(xdot,result, state,[configuration(state);velocity(state)])






mvis = MechanismVisualizer(rob, URDFVisuals(urdf))
open(mvis)
final_time = 100.0

function simple_control!(torques::AbstractVector, t, state::MechanismState)
    q = configuration(state)
    v = velocity(state)


    k = .1
    c = .1

    # c = .1
    # torques[7:9] .= -k*q[8:10] -c*v[7:9]
    torques[7:9] .= -k*q[8:10] - c*v[7:9]


end



closed_loop_dynamics = Dynamics(rob,simple_control!)
controlproblem = ODEProblem(closed_loop_dynamics, state, (0., final_time))

t1 = time()
# sol = DifferentialEquations.solve(controlproblem,RK4(),dt = 1e-4,adaptive = false,progress = true)
sol = DifferentialEquations.solve(controlproblem, Tsit5(),reltol = 1e-10)
# sol = solve(controlproblem, TRBDF2(),reltol = 1e-12)
t2 = time()
delta_t = t2-t1
@show delta_t

# sol = solve(controlproblem, Tsit5(), abs_tol = 1e-10, dt = 0.05)
setanimation!(mvis, sol; realtime_rate = 15)

final_time = 250.0
ts, qs, vs = simulate(state, final_time,simple_control!; Δt = 1e-3,stabilization_gains = nothing)
# ts, qs, vs = simulate(state, final_time)
# sleep(3)
MeshCatMechanisms.animate(mvis, ts, qs; realtimerate = 15.)
