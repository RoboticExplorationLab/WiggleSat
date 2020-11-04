# using Pkg; Pkg.activate(@__DIR__)
using LinearAlgebra
using Attitude
using SatelliteDynamics
using Infiltrator
using FFTW
using MATLAB


function run_SD_prop(dt)


    # Declare simulation initial Epoch
    epc0 = Epoch(2017, 12, 20, 12, 0, 0, 0.0)

    # Declare initial state in terms of osculating orbital elements
    oe0  = [R_EARTH + 550e3, 0.04, 40.0, 0, 0, 0]

    # Convert osculating elements to Cartesean state
    eci0 = sOSCtoCART(oe0, use_degrees=true)

    # Set the propagation end time to one orbit period after the start
    T    = 20*orbit_period(oe0[1])
    epcf = epc0 + T

    # Create an EarthInertialState orbit propagagator
    orb  = EarthInertialState(epc0, eci0, dt=dt,
                area_drag = 0.01, coef_drag = 2.0, area_srp = 0.01, coef_srp = 2.0,
                mass=1.0, n_grav=1, m_grav=1,
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


function eclipse_check(x::Array{<:Real, 1}, r_sun::Array{<:Real, 1})
    """
    Computes the illumination fraction of a satellite in Earth orbit using a
    cylindrical Earth shadow model.

    Arguments:
    - `x::Array{<:Real, 1}`: Satellite Cartesean state in the inertial
                             reference frame [m; m/s]
    - `r_sun::Array{<:Real, 1}`: Position of sun in inertial frame.

    Return:
    - `nu::Float64`: Illumination fraction (0 <= nu <= 1). nu = 0 means
                     spacecraft in complete shadow, nu = 1 mean spacecraft
                     fully illuminated by sun.
    References:
    1. O. Montenbruck, and E. Gill, Satellite Orbits: Models, Methods
                                    and Applications_, 2012, p.80-83.
    """
    # Satellite position ECI
    r = x[1:3]

    # Sun-direction unit-vector
    e_sun = r_sun / norm(r_sun)

    # Projection of spacecraft position
    s = dot(r, e_sun)

    # Compute illumination
    nu = true
    if s/norm(s) >= 1.0 || norm(r - s*e_sun) > R_EARTH
        nu = false
    end

    return nu
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

function srp_torque(r_eci,epoch,ᴺqᴮ,params)

    R_spec = params.R_spec
    R_diff = params.R_diff
    R_abs = params.R_abs

    N_b = params.faces_n
    S = params.faces_areas

    r_sun_eci = sun_position(epoch)

    # position vector from spacecraft to sun
    sc_r_sun = r_sun_eci - r_eci

    # normalize and express in the body frame
    ᴮQᴺ = transpose(dcm_from_q(ᴺqᴮ))
    s = ᴮQᴺ*normalize(sc_r_sun)

    # Get dot products
    I_vec = params.faces_n_b'*s

    F_srp = zeros(6)
    τ_srp = zeros(3)
    for i = 1:6
        F_srp[i] = -P_SUN*S[i]*( 2*(R_diff/3 + R_spec*I_vec[i])*N_b[:,i] +
                    (1-R_spec)*s)*max(I_vec[i],0)
end

# function fft_analysis()
#     x = τ_plot[2,:]
#     y = fft(x)
#     f = (0:length(y)-1)*(1/dt)/length(y)
#     Y_mag = abs.(y)
#     Y_phase = imag(y)
#
#     Y_first_half =
function run_sim()



    faces_n = Array([I(3) -I(3)])
    faces_areas = .01*ones(6)
    R_spec = 0.33
    R_diff = 0.34
    R_abs  = 0.33


    params  = (J = diagm([1;2;3.0]),face_n = faces_n,
               faces_areas = faces_areas, R_spec = R_spec,
               R_diff = R_diff, R_abs = R_abs)
    dt = 20.0

    t, epc, r_eci, v_eci = run_SD_prop(dt)

    ᴺqᴮ = randq()

    N = length(r_eci)

    τ_hist = fill(zeros(3),N)

    for k = 1:N

        τ_hist[k] = gg_torque(r_eci[k],ᴺqᴮ,params)

    end


    τ_plot = mat_from_vec(τ_hist)

    # mat"
    # figure
    # hold on
    # plot($τ_plot')
    # hold off
    # "
    #
    # r_eci = mat_from_vec(r_eci)
    # mat"
    # figure
    # hold on
    # plot3($r_eci(1,:),$r_eci(2,:),$r_eci(3,:) )
    # hold off
    # "
    #
    # mat"
    # figure
    # hold on
    # plot($r_eci')
    # hold off
    # "

    x = τ_plot[2,:]
    y = fft(x)
    f = (0:length(y)-1)*(1/dt)/length(y)
    Y_mag = abs.(y)
    Y_phase = imag(y)

    mat"
    figure
    hold on
    plot($f,$Y_mag)
    %xlim([0 1e-3])
    xlim([0 $(1/dt)/2])
    hold off
    "

    mat"
    figure
    hold on
    plot($τ_plot(2,:))
    hold off"

    mat"
    figure
    hold on
    plot($f,$Y_phase)
    hold off
    "

end

run_sim()
