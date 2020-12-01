using LinearAlgebra


using RigidBodyDynamics
using MeshCat, MeshCatMechanisms

urdf = "/Users/kevintracy/devel/WiggleSat/urdf_multibody/3armsc_6u.xml"


rob = parse_urdf(urdf,gravity = [0.0,0.0,0.0],floating = true)

state = MechanismState(rob)


mvis = MechanismVisualizer(rob, URDFVisuals(urdf))
open(mvis)
final_time = 50.0

ts = [(i-1)*10.0 for i = 1:100]
qs2 = fill(zeros(10),length(ts))
for ii =1:length(ts)
    qs2[ii] = [1;0;0;0;0;0;0;0]
end

MeshCatMechanisms.animate(mvis, ts, qs2; realtimerate = 500.)
