using Pkg, LinearAlgebra

cd("/Users/kevintracy/julia_research/urdf_multibody/RBD_env")
Pkg.activate(".")

using RigidBodyDynamics, RigidBodySim
using MeshCat, MeshCatMechanisms

urdf = "/Users/kevintracy/julia_research/urdf_multibody/3armsc.xml"

# rob = parse_urdf(urdf,floating = true)
rob = parse_urdf(urdf,gravity = [0.0,0.0,0.0],floating = true)

state = MechanismState(rob)

theta0 = deg2rad(30)
# r0 = normalize([0;0;1.0])
r0 = normalize(randn(3))
q0 = [cos(theta0/2);r0*sin(theta0/2)]
set_configuration!(state,[q0;0;0;0;0;0])
set_velocity!(state,[0;0;0;0;.0;.0;0.0])
# set_velocity!(state,[.05])


# result = DynamicsResult(rob)
#
# dynamics!(result, state,[.1])

mvis = MechanismVisualizer(rob, URDFVisuals(urdf))
open(mvis)
final_time = 50.0

function simple_control!(torques::AbstractVector, t, state::MechanismState)
    q = configuration(state)
    v = velocity(state)

    quat = q[1:4]
    if quat[1]<0.0
        quat .= -quat
    end
    error_phi = quat[2:4];


    # k = .1
    # c = .1
    # torques[7:9] .= -k*q[8:10] -c*v[7:9]

    k = .50
    c = 1.0
    # torques[1:3] .= -k*error_phi - c*v[1:3]
    torques[7:9] .= -(-k*error_phi - c*v[1:3])

    # @show q
    # @infiltrate
    # error()
end

ts, qs, vs = simulate(state, final_time,simple_control!; Δt = 1e-3,stabilization_gains = nothing)

# closed_loop_dynamics = Dynamics(rob,simple_control!)
# controlproblem = ODEProblem(closed_loop_dynamics, state, (0., final_time))

# t1 = time()
# # sol = solve(controlproblem, Vern7(), abs_tol = 1e-7, dt = 0.01)
# sol = solve(controlproblem, Tsit5(), abs_tol = 1e-10)
# t2 = time()
# delta_t = t2-t1
# @show delta_t

# sol = solve(controlproblem, Tsit5(), abs_tol = 1e-10, dt = 0.05)
# setanimation!(mvis, sol; realtime_rate = 5)
# ts, qs, vs = simulate(state, final_time;Δt = 1e-5)

MeshCatMechanisms.animate(mvis, ts, qs; realtimerate = 15.)

# function phi_from_q(q)
#     """axis angle from quaternion (scalar last)"""
#
#     v = q[2:4]
#     s = q[1]
#     normv = norm(v)
#
#     if normv == 0.0
#         return zeros(3)
#     else
#         r = v / normv
#         θ = 2 * atan(normv, s)
#         return r * θ
#     end
# end
#
# phi_errors = zeros(length(qs))
#
# for i = 1:length(qs)
#     phi_errors[i] = norm(phi_from_q(qs[i][1:4]))
# end
#
# mat"
# figure
# hold on
# title('Pointing Error')
# plot(linspace(0,100,length($phi_errors)),rad2deg($phi_errors))
# ylabel('Pointing Error (deg)')
# xlabel('Time (s)')"
