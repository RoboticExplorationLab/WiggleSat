cd("/Users/kevintracy/julia_research/urdf_multibody/RBD_env")
Pkg.activate(".")

using RigidBodyDynamics

urdf = "/Users/kevintracy/julia_research/urdf_multibody/3armsc.xml"

# rob = parse_urdf(urdf,floating = true)
rob = parse_urdf(urdf,gravity = [0.0,0.0,0.0],floating = true)

state = MechanismState(rob)





M = mass_matrix(state)

C = dynamics_bias(state)

vdot =
