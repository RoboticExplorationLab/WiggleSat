

function generate_geometry(length,width,height)
    # dim_L = .1
    # dim_W = .2
    # dim_H = .3
    L = length
    W = width
    H = height
    volume = L*W*H
    mass = volume*(1/0.01)

    # inertia
    J = (mass/12)*Array(I(3))
    J[1,1] *= H^2 +W^2
    J[2,2] *= H^2 +L^2
    J[3,3] *= L^2 +W^2

    # faces normal vectors
    faces_n = vec_from_mat(Array([I(3) -I(3)]))

    # faces areas
    faces_areas = [W*H; H*L; L*W; W*H; H*L; L*W]

    # face position vectors
    r = fill(zeros(3),6)
    r[1] = L/2 * faces_n[1]
    r[2] = W/2 * faces_n[2]
    r[3] = H/2 * faces_n[3]
    r[4] = L/2 * faces_n[4]
    r[5] = W/2 * faces_n[5]
    r[6] = H/2 * faces_n[6]

    # properties assuming 1/3 aluminum, 2/3 gallium arsenide
    # R_spec = 0.3
    # R_diff = .23333
    # R_abs = 0.46667
    # R_spec = 0.1
    # R_diff = .2
    # R_abs = 0.7
    R_spec = 0.0
    R_diff = 0.0
    R_abs = 1.0

    faces = (normal_vecs = faces_n, areas = faces_areas, r = r,
             R_spec = R_spec, R_diff = R_diff, R_abs  = R_abs)

    return faces, J
end





faces, J = generate_geometry(.1,.2,1.3)

params  = (J = J,faces = faces)

r_eci = [6.651010847999999e6;0;0]
epoch = Epoch(2017, 12, 20, 12, 0, 0, 0.0)
ᴺqᴮ = [0;0;0;1]

R_spec = params.faces.R_spec
R_diff = params.faces.R_diff
R_abs = params.faces.R_abs

N_b = params.faces.normal_vecs
S = params.faces.areas
r = params.faces.r

r_sun_eci = sun_position(epoch);
# r_sun_eci = norm(r_sun_eci)*normalize([1;1;0])

# position vector from spacecraft to sun
sc_r_sun = (r_sun_eci - r_eci)

# normalize and express in the body frame
ᴮQᴺ = transpose(dcm_from_q(ᴺqᴮ))
s = ᴮQᴺ*normalize(sc_r_sun)
# @show s
# Get dot products
# I_vec = N_b'*s
# @show I_vec

global F_srp = fill(zeros(3),6)
global τ_srp = zeros(3)
for i = 1:6
    # @show i
    # @infiltrate
    # error()
    n = N_b[i]
    # @show norm(s)
    # @show norm(n)

    cosθ = dot(n,s)

    F_srp[i] = -P_SUN*S[i]*( 2*(R_diff/3 + R_spec*cosθ)*n +
                (1-R_spec)*s)*max(cosθ,0.0)

    # @show F_srp
    global τ_srp += cross(r[i],F_srp[i])
    @show τ_srp
    # @infiltrate
end
# error()
# @show τ_srp
active_faces = sum(I_vec .> 0)

@show τ_srp
@show F_srp
