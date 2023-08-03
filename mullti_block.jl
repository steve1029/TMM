include("./TMM_module.jl")
using .TMM
using Base
using Distributions
using LinearAlgebra
using NPZ
using Printf
using Plots

# Choose the backends of the Plots package.
gr()
# pythonplot()
# plotlyjs()
# unicodeplots()

um = 1e-6
nm = 1e-9
λ0 = 300 * nm # wavelength in free space.
k0 = 2*π / λ0 # wavenumber in free space.
c = 299792458 # m/s.
freq = c / λ0 # frequency
ω = 2π * freq
μ_0 = 4*π*10^-7
ε_0 = 8.8541878128e-12

# polar angle and azimuthal angle.
# θ = π / 12 # radian.
# ϕ = π / 20 # radian.
θ = -π / 6
# θ = 0
ϕ = 0
# ϕ = π / 6

dx = 10*nm
dz = 10*nm

Lx = 1000*nm
Lz = 2000*nm

z1 =  500*nm
z2 =  700*nm
z3 = 1000*nm
z4 = 1200*nm
z5 = 1500*nm

N = 4
zs = [z1, z2, z3, z4, z5]
murs = [1, 1, 1, 1]
epsrs = [1.3^2, 1.7^2, 1.5^2, 1.2^2]
imp = sqrt(μ_0 /ε_0)
input = [0., 1.] # [Eix, Eiy]

Sn, Can, Cbn = redheffer_n_blocks(N, μ_0, ω, imp, k0, θ, ϕ, zs, murs, epsrs)

println(size(Sn))
println(size(Can), typeof(Can))
println(size(Cbn), typeof(Cbn))
println(size(Can[1]))
println(size(Can[2]))

"""
# Get the S matrix of a left half-infinite block.
eigvals0, eigvecs0 = TMM.get_eigenvectors(k0, kx0, ky0, mur, epsr, impedance)
Wh, Vh = TMM._make_WhVh(μ_0, ω, kx0, ky0, kz0)
Wp, Vp = TMM._make_WpVp(eigvals0, eigvecs0, zm)
Wm, Vm = TMM._make_WmVm(eigvals0, eigvecs0, zm)

A = inv(Wh)*Wm - inv(Vh)*Vm
B = inv(Wh)*Wp - inv(Vh)*Vp
C = inv(Wm)*Wh - inv(Vm)*Vh
D = inv(Wm)*Wp - inv(Vm)*Vp
E = inv(Wm)*Wh + inv(Vm)*Vh

R00_RtoL =  -inv(A) * B 
T00_LtoR =   inv(C) * D
T00_RtoL = 2*inv(A)
R00_LtoR =  -inv(C) * E

S00 = zeros(ComplexF64, 4, 4)
S00[1:2, 1:2] = T00_LtoR
S00[1:2, 3:4] = R00_LtoR
S00[3:4, 1:2] = R00_RtoL
S00[3:4, 3:4] = T00_RtoL

# Get the S matrix of a right half-infinite block.
"""