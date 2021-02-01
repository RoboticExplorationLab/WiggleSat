using LinearAlgebra
using Statistics
using JLD2
using MATLAB
using FiniteDiff
using StaticArrays
using SparseArrays
using ProgressMeter
@load "/Users/kevintracy/devel/WiggleSat/orbit_data_zac.jld2" τ_hist B_hist_b J


function get_params(τ_hist, J)
    dt = 1.0 #seconds
    #t = Array(LinRange(0.0, dt*length(τ_hist), length(τ_hist)+1))
    t = Array(1:600)

    J_boom = Diagonal(SVector{3}(0.25.*[.1^2; .2^2; .3^2]))
    rad2arcsec = 3600*180/pi

    Nm2μNm = 1e6
    τ = Nm2μNm.*hcat(τ_hist...)

    #Double integrator dynamics
    #Units are arcsec, seconds, torque in μNm
    #State is [θ; θ̇]
    A = [zeros(3,3) I; zeros(3,6)]
    B = Array([zeros(3,3); (rad2arcsec/Nm2μNm).*inv(J)]) #control input jacobian
    # C = Array(Diagonal(ones(6)))
    C = Diagonal(@SVector ones(6))
    G = [zeros(3,3); (rad2arcsec/Nm2μNm).*inv(J)] #disturbance input jacobian

    #Convert to discrete time
    H = exp(dt.*[A B G; zeros(6,12)])
    A = SMatrix{6,6}(H[1:6,1:6])
    B = SMatrix{6,3}(H[1:6,7:9])
    G = SMatrix{6,3}(H[1:6,10:12])

    #Kalman Filter Design
    # V = [0.00001*I zeros(6,3); zeros(3,6) 0.0001*I] # TODO: tune this
    #
    #
    # W_telescope = 1.0^2 #arcsec^2 1-sigma at 1 Hz
    # #W_gyro = (0.06*60)^2 #(arcsec/sec)^2 1-sigma at 1 Hz for Epson MEMS IMU (ARW = 0.06 deg/sqrt(hr))
    # W_gyro = (0.0035*60)^2 #(arcsec/sec)^2 1-sigma at 1 Hz for Honeywell GG1320 laser gyro (ARW = 0.0035 deg/sqrt(hr))
    # W = Array(Diagonal([(W_telescope/3.0)*ones(3); (W_gyro/3.0)*ones(3)])) #

    #Double integrator dynamics + bias torque
    #Units are arcsec, seconds, torque in μNm
    #State is [θ; θ̇; τ_b]
    Af = [zeros(3,3) I zeros(3,3); zeros(3,6) (rad2arcsec/Nm2μNm).*inv(J); zeros(3,9)]
    Bf = Array([zeros(3,3); (rad2arcsec/Nm2μNm).*inv(J); zeros(3,3)]) #control input jacobian
    Cf = sparse([Array(Diagonal(ones(6))) zeros(6,3)])

    #Convert to discrete time
    Hf = exp(dt.*[Af Bf; zeros(3,12)])
    Af = SMatrix{9,9}(Hf[1:9,1:9])
    Bf = SMatrix{9,3}(Hf[1:9,10:12])
    # @infiltrate
    # error()

    return τ, A, B, C, G, Af, Bf, Cf
end
function run_zac(τ_hist,B_hist_b,J,τ, A, B, C, G, Af, Bf, Cf,W_telescope ,β)
dt = 1.0 #seconds
#t = Array(LinRange(0.0, dt*length(τ_hist), length(τ_hist)+1))
t = Array(1:600)

J_boom = Diagonal(SVector{3}(0.25.*[.1^2; .2^2; .3^2]))
rad2arcsec = 3600*180/pi

Nm2μNm = 1e6
τ = Nm2μNm.*hcat(τ_hist...)

#Double integrator dynamics
#Units are arcsec, seconds, torque in μNm
#State is [θ; θ̇]
# A = [zeros(3,3) I; zeros(3,6)]
# B = Array([zeros(3,3); (rad2arcsec/Nm2μNm).*inv(J)]) #control input jacobian
# C = Array(Diagonal(ones(6)))
# G = [zeros(3,3); (rad2arcsec/Nm2μNm).*inv(J)] #disturbance input jacobian
#
# #Convert to discrete time
# H = exp(dt.*[A B G; zeros(6,12)])
# A = H[1:6,1:6]
# B = H[1:6,7:9]
# G = H[1:6,10:12];
# τ, A, B, C, G, Af, Bf, Cf =
#Kalman Filter Design
# V = [0.00001*I zeros(6,3); zeros(3,6) 0.0001*I] # TODO: tune this
# V = [β[1]*I zeros(6,3); zeros(3,6) β[2]*I] # TODO: tune this
V = Diagonal( SA[β[1],β[1],β[1],β[6],β[6],β[6],β[2],β[2],β[2]])

# W_telescope = 1.0^2 #arcsec^2 1-sigma at 1 Hz # NOTE: THIS IS NOW AN INPUT
#W_gyro = (0.06*60)^2 #(arcsec/sec)^2 1-sigma at 1 Hz for Epson MEMS IMU (ARW = 0.06 deg/sqrt(hr))
W_gyro = (0.0035*60)^2 #(arcsec/sec)^2 1-sigma at 1 Hz for Honeywell GG1320 laser gyro (ARW = 0.0035 deg/sqrt(hr))
# W = Array(Diagonal([(W_telescope/3.0)*ones(3); (W_gyro/3.0)*ones(3)])) #
W = Diagonal(SA[W_telescope/3,W_telescope/3,W_telescope/3,(W_gyro/3.0),(W_gyro/3.0),(W_gyro/3.0)])
#Double integrator dynamics + bias torque
#Units are arcsec, seconds, torque in μNm
#State is [θ; θ̇; τ_b]
# Af = [zeros(3,3) I zeros(3,3); zeros(3,6) (rad2arcsec/Nm2μNm).*inv(J); zeros(3,9)]
# Bf = Array([zeros(3,3); (rad2arcsec/Nm2μNm).*inv(J); zeros(3,3)]) #control input jacobian
# Cf = [Array(Diagonal(ones(6))) zeros(6,3)]
#
# #Convert to discrete time
# Hf = exp(dt.*[Af Bf; zeros(3,12)])
# Af = Hf[1:9,1:9]
# Bf = Hf[1:9,10:12]

#LQR Controller Design with Bias Torque
#Note that dlqr doesn't work since the system is not strictly controllable
# Q = Array(Diagonal([1.0*ones(3); 10.0*ones(3); zeros(3)])) # NOTE: tune this
# R = Array(Diagonal(0.01.*ones(3)))                         # NOTE: tune this
Q = (Diagonal([β[3]*ones(3); β[4]*ones(3); zeros(3)])) # NOTE: tune this
R = (Diagonal(β[5].*ones(3)))

P = copy(Array(Q))
K = zeros(3,9)
for k = 1:100
    K = (R+Bf'*P*Bf)\Bf'*P*Af
    P = Q + K'*R*K + (Af-Bf*K)'*P*(Af-Bf*K)
end

#Closed-loop sim
x = zeros(6,length(t))
u = zeros(3,length(t)-1)
# initial condition
# x[:,1] .= [1.0*randn(3); 0.1*randn(3)]
# x[:,1] .= [0*randn(3); 0*randn(3)]

x̄ = zeros(9,length(t)) #Filter state includes torque bias
x̄[:,1] = [zeros(6);τ[:,1]]
P = zeros(9,9,length(t))
P[:,:,1] .= Array(Diagonal([100*ones(6); ones(3)]))

sqrtW = sqrt(W)
L = zeros(size(x̄,1),6)
for k = 1:(length(t)-1)
    #Run controller one step
    u[:,k] .= -K*x̄[:,k] #[x[:,k]; τ[:,k]]

    # true dynamics update
    x[:,k+1] .= A*x[:,k] + B*u[:,k] + G*τ[:,k]

    # true measurement update
    y = C*x[:,k+1] + sqrtW*randn(6) #Generate measurement with appropriate noise

    # KF predict
    xp = Af*x̄[:,k] + Bf*u[:,k] #State prediction
    Pp = Af*P[:,:,k]*Af' + V #Prediction covariance

    # KF innovate
    z = y - Cf*xp #Innovation
    S = Cf*Pp*Cf' + W #Innovation covariance
    # L = Pp*Cf'*inv(S) #Kalman gain
    # in case the weights cause this to fail
    try
        L = Pp*Cf'/S #Kalman gain
    catch
        L *=NaN
    end

    # KF update
    x̄[:,k+1] .= xp + L*z #Measurement update
    P[:,:,k+1] .= Pp - L*S*L' #Covariance update
    for m = 1:9
        P[m,m,k+1] = max(1e-6,P[m,m,k+1])
    end
end


Nt = length(t)
Xvar = 0.0
for k = 1:Nt
    Xvar += (1/Nt)*(x[1:3,k]'*x[1:3,k])
end
RMSerr = sqrt(Xvar)
return RMSerr
# @show RMSerr
# n = 2
# plot(x[n,:])
# plot!(x̄[n,:])

# n = 7
# plot(x̄[n,t])
# plot!(τ[n-6,t])

# # plot(x[1,:],x[2,:])
# xcirc = zeros(361)
# ycirc = zeros(361)
# for k = 0:360
#     xcirc[k+1] = cosd(k)
#     ycirc[k+1] = sind(k)
# end
# # plot!(xcirc,ycirc,linewidth=2)
# # xlims!(-2,2)
# # ylims!(-2,2)
#
#
# mat"""
# plot($x(1,:),$x(2,:),'linewidth',2)
# hold on
# plot($xcirc,$ycirc,'linewidth',2)
# xlim([-1.5 1.5])
# ylim([-1.5 1.5])
# axis equal
# xlabel('Yaw (arcsec)')
# ylabel('Pitch (arcsec)')
# %addpath('~/Documents/MATLAB/matlab2tikz/src')
# %matlab2tikz('closed-loop-pointing.tikz')
# """
#
# 0.06*60

end


function driver(τ_hist,B_hist_b,J,W_telescope ,β)
τ, A, B, C, G, Af, Bf, Cf = get_params(τ_hist, J)
trials = 20
rms_vec = zeros(trials)
for i = 1:trials
rms_vec[i] = run_zac(τ_hist,B_hist_b,J,τ, A, B, C, G, Af, Bf, Cf,W_telescope ,β)
end
return mean(rms_vec)
end

function optimize_pointing(τ_hist,B_hist_b,J,β_0,W_telescope,N,α)
# β_0 = [0.00001,0.0001,1,10,0.01]
# β_0 = [1.1671759312723126e-8, 7.190344939687723e-5, 1.6697629325758443, 0.057106137532170496, 0.0002467539801031718]

# W_telescope = 1.0^2 #arcsec^2 1-sigma at 1 Hz
fd_closure(βc) = driver(τ_hist,B_hist_b,J,W_telescope,βc)
bestβ = β_0
bestJ = fd_closure(bestβ)
@showprogress "optimizing..." for i = 1:N
    β = (abs.(1 .+ α*randn(6))) .* bestβ
    # β = copy(bestβ)
    # β[rand(1:6)] *= (abs(1 + α*randn()))
    newJ =  fd_closure(β)
    println(bestJ)
    if newJ < bestJ
        bestJ = (newJ)
        bestβ = copy(β)
    end
end

return bestβ, bestJ
end

# # 1 arcsecond case
# W_telescope = 1^2 #arcsec^2 1-sigma at 1 Hz
# β_0 = [7.858862760110066e-7, 6.744841534181973e-5, 91.72214577240014,
# 0.19447062964953343, 0.007291584771687058, 7.287578908992002e-7]
# # β_0 = [0.00001,0.0001,1,10,0.01,0.00001]
# bestβ_1,bestJ_1 = optimize_pointing(τ_hist,B_hist_b,J,β_0,W_telescope,100,1)
#
# # .1 arcsecond case
# W_telescope = 0.1^2 #arcsec^2 1-sigma at 1 Hz
# β_0 = [4.332303216588739e-7, 6.990003715236112e-5, 42.12592959061058,
#  0.1404486745411262, 0.0002723656503952372, 1.6512704231365986e-7]
# # β_0 = [0.00001,0.0001,1,10,0.01,0.00001]
# bestβ_1,bestJ_1 = optimize_pointing(τ_hist,B_hist_b,J,β_0,W_telescope,100,1)

# # .5 arcsecond case
# W_telescope = 0.5^2 #arcsec^2 1-sigma at 1 Hz
# β_0 = [7.586731861311927e-9, 0.0005854287799096624, 1.0881005224942581,
#  0.018258319791698032, 0.0007999924517238642, 1.4672777505792572e-9]
# # β_0 = [0.00001,0.0001,1,10,0.01,0.00001]
# bestβ_1,bestJ_1 = optimize_pointing(τ_hist,B_hist_b,J,β_0,W_telescope,100,1)

# .1 arcsecond case
W_telescope = 0.05^2 #arcsec^2 1-sigma at 1 Hz
β_0 = [4.332303216588739e-7, 6.990003715236112e-5, 42.12592959061058,
 0.1404486745411262, 0.0002723656503952372, 1.6512704231365986e-7]
# β_0 = [0.00001,0.0001,1,10,0.01,0.00001]
bestβ_1,bestJ_1 = optimize_pointing(τ_hist,B_hist_b,J,β_0,W_telescope,100,1)







# # 10 arcsecond case
# W_telescope = 10.0^2 #arcsec^2 1-sigma at 1 Hz
# β_0 = [8.567935963068527e-11, 2.4282618041065227e-5,
# 0.03512714908159575, 0.13183770364640643, 0.00021472969710868088]
# bestβ_2,bestJ_2 = optimize_pointing(τ_hist,B_hist_b,J,β_0,W_telescope,100,1)

# # 10 arcsecond case
# W_telescope = 0.1^2 #arcsec^2 1-sigma at 1 Hz
# β_0 = [4.859209940424407e-10, 0.00017843343145081865, 2.004553423591832,
#  2.3027986802413978e-5, 5.574007444038023e-5]
# bestβ_3,bestJ_3 = optimize_pointing(τ_hist,B_hist_b,J,β_0,W_telescope,100,.001)


arc_sec_vec = [0.1,  0.5,  1]
rms_vec     = [0.33, 0.39, .49 ]
