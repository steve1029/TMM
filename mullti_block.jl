include("./TMM_module.jl")
using .single_block
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
λ = 300 * nm
k0 = 2*π / λ # wavenumber in free space.
c = 299792458 # m/s.
μ_0 = 4*π*10^-7
ε_0 = 8.8541878128e-12
mur = 1
epsr = 2.0^2
impedance = sqrt(μ_0 /ε_0)

# polar angle and azimuthal angle.
# θ = π / 12 # radian.
# ϕ = π / 20 # radian.
θ = -π / 6
# θ = 0
ϕ = 0
# ϕ = π / 6

# wave vector in free space.
kx0 = k0 * sin(θ) * cos(ϕ)
ky0 = k0 * sin(θ) * sin(ϕ)
kz0 = k0 * cos(θ)

dx = 10*nm
dz = 10*nm

Lx = 1000*nm
Lz = 1500*nm

zm =  600*nm
zp = 1000*nm

input = [0., 1.] # [Eix, Eiy]

# Get the S matrix of a left half-infinite block.
eigvals0, eigvecs0 = get_eigenvectors(k0, kx0, ky0, mur, epsr, impedance)
Wh, Vh = make_WhVh(μ_0, ω, kx0, ky0, kz0)
Wp, Vp = make_WpVp(eigvals, eigvecs, zm)

S = zeros(ComplexF64, 4, 4)

S[1:2, 1:2] = aT
S[1:2, 3:4] = aR
S[3:4, 1:2] = bR
S[3:4, 3:4] = bT
