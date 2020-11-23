using Pkg, LinearAlgebra

cd(joinpath(@__DIR__,"RBD_env"))
Pkg.activate(".")

using RigidBodyDynamics
using MeshCat, MeshCatMechanisms
using JLD2
@load "/Users/kevintracy/devel/WiggleSat/state_hist_for_vis.jld2" θm
# urdf = "/Users/kevintracy/julia_research/urdf_multibody/3armsc.xml"
urdf = joinpath(@__DIR__,"3armsc.xml")


rob = parse_urdf(urdf,gravity = [0.0,0.0,0.0],floating = true)

mvis = MechanismVisualizer(rob, URDFVisuals(urdf))
open(mvis)

# ts = 0:10.0:((length(θm)-1)*10)
ts = [(i-1)*10.0 for i = 1:length(θm)]
qs2 = fill(zeros(10),length(ts))
for ii =1:length(ts)
    qs2[ii] = [1;0;0;0;0;0;0;θm[ii]]
    # qs2[ii][1:4] = [1.0; 0; 0; 0]
    # qs2[ii][8:10] = θm[ii]
    # @show θm[ii]
end

MeshCatMechanisms.animate(mvis, ts, qs2; realtimerate = 500.)



X = fill(zeros(4),10)
for i = 1:10
    X[i][1:2] = [3;4.0]
    X[i][3:4] = randn(2)
end
