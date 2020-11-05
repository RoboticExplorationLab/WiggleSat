# using Pkg; Pkg.activate(@__DIR__)
using LinearAlgebra
using Attitude
using SatelliteDynamics
using Infiltrator
using FFTW
using MATLAB


function run_SD_prop(dt)


    # Declare simulation initial Epoch
    epc0 = Epoch(2020, 5, 1, 1, 0, 0, 0.0)

    # Declare initial state in terms of osculating orbital elements
    oe0  = [R_EARTH + rand_in_range(300.0,800.0)*1e3, rand_in_range(0,0.01),
            rand_in_range(0.0,89.9), 0, 0, 0]

    # Convert osculating elements to Cartesean state
    eci0 = sOSCtoCART(oe0, use_degrees=true)

    # Set the propagation end time to one orbit period after the start
    T    = 0.2*orbit_period(oe0[1])
    epcf = epc0 + T

    # Create an EarthInertialState orbit propagagator
    orb  = EarthInertialState(epc0, eci0, dt=dt,
                area_drag = 0.01, coef_drag = 2.0, area_srp = 0.01, coef_srp = 2.0,
                mass=1.0, n_grav=2, m_grav=2,
                drag=false, srp=false,
                moon=false, sun=false,
                relativity=false
    )

    # Propagate the orbit
    t, epc, eci = sim!(orb, epcf)

    r_eci = vec_from_mat(eci[1:3,:])
    v_eci = vec_from_mat(eci[4:6,:])
    return t, epc, r_eci, v_eci
end


function gg_torque(r_eci,ᴺqᴮ,params)
    """Environmental torque modeling.

    Args:


    Output: torque (N⋅m)

    """

    # distance from center of earth
    r = norm(r_eci)

    # normalized position vector expressed in body basis
    ᴮQᴺ = transpose(dcm_from_q(ᴺqᴮ))
    n = normalize(ᴮQᴺ*r_eci)

    τ_gg = (3*GM_EARTH/(r^3))*cross(n,params.J*n)

    return τ_gg
end

function aerodynamic_torque(r_eci,v_eci,epoch,ᴺqᴮ,params)

    Cd = params.faces.Cd
    S = params.faces.areas
    r = params.faces.r

    r_sun_eci = sun_position(epoch)
    ρ = density_harris_priester([r_eci;v_eci],r_sun_eci)

    ᴮQᴺ = transpose(dcm_from_q(ᴺqᴮ))

    v_rel_eci = v_eci - cross([0;0;OMEGA_EARTH],r_eci)

    v_rel_b = ᴮQᴺ*v_rel_eci

    v = normalize(v_rel_b)

    F_aero = fill(zeros(3),6)
    τ_aero = zeros(3)

    for i = 1:6
        n = params.faces.normal_vecs[:,i]

        cosθ = dot(n,v)

        F_aero[i] = -.5*ρ*Cd*norm(v_rel_b)*v_rel_b*S[i]*max(cosθ,0)

        τ_aero += cross(r[i],F_aero[i])
    end

    return τ_aero
end


function srp_torque(r_eci,epoch,ᴺqᴮ,params)

    R_spec = params.faces.R_spec
    R_diff = params.faces.R_diff
    R_abs = params.faces.R_abs

    N_b = params.faces.normal_vecs
    S = params.faces.areas
    r = params.faces.r

    r_sun_eci = sun_position(epoch)

    # position vector from spacecraft to sun
    sc_r_sun = (r_sun_eci - r_eci)

    # normalize and express in the body frame
    ᴮQᴺ = transpose(dcm_from_q(ᴺqᴮ))
    s = ᴮQᴺ*normalize(sc_r_sun)


    F_srp = fill(zeros(3),6)
    τ_srp = zeros(3)
    for i = 1:6

        n = N_b[:,i]
        cosθ = dot(n,s)

        F_srp[i] = -P_SUN*S[i]*( 2*(R_diff/3 + R_spec*cosθ)*n +
                    (1-R_spec)*s)*max(cosθ,0.0)

        τ_srp += cross(r[i],F_srp[i])

    end


    return eclipse_conical(r_eci,r_sun_eci)*τ_srp
    # return τ_srp
end

function sample_inertia(J,deg_std,percent_std)
    E = eigen(J)
    S = @views E.vectors
    D = @views E.values
    # scale moments
    # @infiltrate
    # error()
    J_new = S * Diagonal(((1 .+(percent_std/100)*randn(3)) .* D)) * S'

    # J_new = E.vectors*
    #         Diagonal(((1 .+(percent_std/100)*randn(3)) .* E.values))
    #         *E.vectors'

    # rotate
    R = dcm_from_phi(deg2rad(deg_std)*randn(3))
    J_new = R*J_new*R'

    return J_new
end




function generate_geometry(length,width,height)

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
    J = sample_inertia(J,20.0,20.0)

    # faces normal vectors
    faces_n = Array([I(3) -I(3)])

    # faces areas
    faces_areas = [W*H; H*L; L*W; W*H; H*L; L*W]

    # CG offset
    scale_factor = 2
    cg_offset = [rand_in_range(-L/scale_factor,L/scale_factor);
                 rand_in_range(-W/scale_factor,W/scale_factor);
                 rand_in_range(-H/scale_factor,H/scale_factor)]
    # face position vectors
    r = fill(zeros(3),6)
    r[1] = L/2 * faces_n[:,1] - cg_offset
    r[2] = W/2 * faces_n[:,2] - cg_offset
    r[3] = H/2 * faces_n[:,3] - cg_offset
    r[4] = L/2 * faces_n[:,4] - cg_offset
    r[5] = W/2 * faces_n[:,5] - cg_offset
    r[6] = H/2 * faces_n[:,6] - cg_offset

    # properties assuming 1/3 aluminum, 2/3 gallium arsenide
    R_spec = 0.3
    R_diff = .23333
    R_abs = 0.46667

    faces = (normal_vecs = faces_n, areas = faces_areas, r = r,
             R_spec = R_spec, R_diff = R_diff, R_abs  = R_abs,
             Cd = 1.5)

    return faces, J
end


function create_spline(t,x)

    spl = Spline1D(t,x)
    # @infiltrate
    # error()
    newt = 0:1:t[end]
    newy = spl(newt)
    return newy
end


function fft_analysis(x,dt)

    # FFTW dft
    y = fft(x)

    # frequencies
    f = (0:length(y)-1)*(1/dt)/length(y)

    # get magnitude at each frequency
    Y_mag = abs.(y)

    return Y_mag, f
end

function run_sim()


    # create spacecraft geometrical properties
    faces, J = generate_geometry(.3,.1,.1)

    params  = (J = J,faces = faces)
    dt = 30.0

    t, epc, r_eci, v_eci = run_SD_prop(dt)

    ᴺqᴮ = randq()
    # ᴺqᴮ = [0;0;0;1]

    N = length(r_eci)

    τ_hist = fill(zeros(3),N)
    B_hist_b = fill(zeros(3),N)
    F_srp = fill(  fill(zeros(3),6)   ,N)
    active_faces = zeros(N)

    for k = 1:N
        τ_hist[k] += gg_torque(r_eci[k],ᴺqᴮ,params)
        τ_hist[k] += srp_torque(r_eci[k],epc[k],ᴺqᴮ,params)
        τ_hist[k] += aerodynamic_torque(r_eci[k],v_eci[k],epc[k],ᴺqᴮ,params)

        B_hist_b[k] = dcm_from_q(ᴺqᴮ)'*IGRF13(r_eci[k],epc[k])
    end


    τ_plot = mat_from_vec(τ_hist)
    # F_plot = mat_from_vec(F_srp)

    τ1 = create_spline(t,τ_plot[1,:])
    τ2 = create_spline(t,τ_plot[2,:])
    τ3 = create_spline(t,τ_plot[3,:])
    dt = 1.0
    Y_mag3, f = fft_analysis(τ1,dt)
    Y_mag2, f = fft_analysis(τ2,dt)
    Y_mag1, f = fft_analysis(τ3,dt)
    Y_mag = Y_mag1 + Y_mag2 + Y_mag3
    # x = τ_plot[2,:]
    # y = fft(x)
    # f = (0:length(y)-1)*(1/dt)/length(y)
    # Y_mag = abs.(y)
    # Y_phase = imag(y)

    # mat"
    # figure
    # hold on
    # plot($f(2:end),$Y_mag(2:end))
    # %xlim([0 1e-3])
    # %xlim([0 $(1/dt)/2])
    # xlim([0 5e-3])
    # ylabel('Magnitude')
    # xlabel('Frequency (Hz)')
    # hold off
    # "
    #
    # mat"
    # figure
    # hold on
    # plot($τ_plot(2,:))
    # hold off"

    # mat"
    # figure
    # hold on
    # plot($f,$Y_phase)
    # hold off
    # "
    return τ_hist, B_hist_b

end

τ_hist, B_hist_b  = run_sim()


# using JuMP
# using OSQP
# N = length(τ_hist)
#
# model = Model(OSQP.Optimizer)
#
# @variable(model, τ[1:3,1:N])
# @variable(model, m[1:3,1:N])
#
# for i = 1:N
#     @constraint(model,con, skew_from_vec(B_hist_b[i])*m[:,i] + τ[:,i] .== τ_hist[i])
# end

using Convex, COSMO

τ = Variable(3,N)
m = Variable(3,N)

cons = Constraint


# cons = Constraint[ ]
for i = 1:N
    push!(cons, skew_from_vec(B_hist_b[i])*m[:,i] + τ[:,i] == τ_hist[i]  )
end
