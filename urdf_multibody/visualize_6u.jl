using LinearAlgebra

using CoordinateTransformations
using RigidBodyDynamics
using MeshCat, MeshCatMechanisms
using GeometryTypes:
    GeometryTypes, HyperRectangle, Vec, Point,
    HomogenousMesh, SignedDistanceField, HyperSphere, GLUVMesh, Pyramid

urdf = "/Users/kevintracy/devel/WiggleSat/urdf_multibody/3armsc_6u.xml"


rob = parse_urdf(urdf,gravity = [0.0,0.0,0.0],floating = true)

state = MechanismState(rob)


mvis = MechanismVisualizer(rob, URDFVisuals(urdf))
# setprop!(mvis["/Background"], "top_color", colorant="red")
# delete!(mvis["/Background"])
setprop!(mvis["/Lights/AmbientLight"],"intensity",10.0)
open(mvis)


final_time = 50.0

r = [0;0;1]
θ = deg2rad(180)
# ϕ = r*θ
quat = [cos(θ/2);r*sin(θ/2)]
ts = [(i-1)*10.0 for i = 1:100]
qs2 = fill(zeros(10),length(ts))
for ii =1:length(ts)
    qs2[ii] = [quat;0;0;0;-1.4;2.5;0.94]
end


MeshCatMechanisms.animate(mvis, ts, qs2; realtimerate = 500.)

cubesat_dims = [3.1;1.01;1.1]
target_cubesat_image = PngImage("/Users/kevintracy/devel/WiggleSat/urdf_multibody/white_solar_panel.png")
target_cubesat_texture = Texture(image=target_cubesat_image)
target_cubesat_material = MeshLambertMaterial(map=target_cubesat_texture)
target_cubesat = HyperRectangle(-Vec(cubesat_dims...)./2, Vec(cubesat_dims...))
setobject!(mvis["target_cubesat"], target_cubesat, target_cubesat_material)
R = Array(float(I(3)))
target_rotation = LinearMap(RotMatrix(R...))
d = [0;0.5;0]
target_translation = Translation(d...)
settransform!(mvis["target_cubesat"],compose(target_translation,target_rotation))


cubesat_dims = [3.1;1.01;1.1]
target_cubesat_image = PngImage("/Users/kevintracy/devel/WiggleSat/urdf_multibody/white_solar_panel.png")
target_cubesat_texture = Texture(image=target_cubesat_image)
target_cubesat_material = MeshLambertMaterial(map=target_cubesat_texture)
target_cubesat = HyperRectangle(-Vec(cubesat_dims...)./2, Vec(cubesat_dims...))
setobject!(mvis["target_cubesat2"], target_cubesat, target_cubesat_material)
R = Array(float(I(3)))
target_rotation = LinearMap(RotMatrix(R...))
# target_rotation = LinearMap(RotZ(θ))
d = [0;-0.5;0]
target_translation = Translation(d...)
settransform!(mvis["target_cubesat2"],compose(target_translation,target_rotation))

cubesat_dims = [0.1;1.0;2.0]
target_cubesat_image = PngImage("/Users/kevintracy/devel/WiggleSat/urdf_multibody/front.png")
target_cubesat_texture = Texture(image=target_cubesat_image)
target_cubesat_material = MeshLambertMaterial(map=target_cubesat_texture)
target_cubesat = HyperRectangle(-Vec(cubesat_dims...)./2, Vec(cubesat_dims...))
setobject!(mvis["target_cubesat3"], target_cubesat, target_cubesat_material)
# R = Array(float(I(3)))
θ = pi/2
target_rotation = LinearMap(RotX(θ))
d = [-1.51;0;0]
target_translation = Translation(d...)
settransform!(mvis["target_cubesat3"],compose(target_translation,target_rotation))
