using LinearAlgebra
using RigidBodyDynamics, RigidBodySim
using MeshCat, MeshCatMechanisms

# path to URDF, yours will be different
#TODO: you guys have to update this path to be wherever you put the urdf file
urdf = "/Users/kevintracy/julia_research/urdf_multibody/3armsc.xml"

# create a robot out of this urdf
rob = parse_urdf(urdf,gravity = [0.0,0.0,0.0],floating = true)

# turn it into a mechanism
state = MechanismState(rob)

# set initial orientation with a quaternion
theta0 = deg2rad(30)
r0 = normalize(randn(3))
q0 = [cos(theta0/2);r0*sin(theta0/2)]

# set the state
set_configuration!(state,[q0;0;0;0;0;0])
set_velocity!(state,[0;0;0;0;.0;.0;0.0])

# start visualization
mvis = MechanismVisualizer(rob, URDFVisuals(urdf))
open(mvis)

# simulate for 50 seconds
final_time = 50.0

function simple_control!(torques::AbstractVector, t, state::MechanismState)
    """This is a simple PD controller for controlling the attitude"""
    q = configuration(state)
    v = velocity(state)

    quat = q[1:4]
    if quat[1]<0.0
        quat .= -quat
    end
    error_phi = quat[2:4];

    k = .50
    c = 1.0

    torques[7:9] .= -(-k*error_phi - c*v[1:3])

end

ts, qs, vs = simulate(state, final_time,simple_control!; Δt = 1e-3,stabilization_gains = nothing)

MeshCatMechanisms.animate(mvis, ts, qs; realtimerate = 15.)
