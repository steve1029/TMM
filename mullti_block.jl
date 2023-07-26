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
mur = 1
epsr = 2.0^2
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
Lz = 1500*nm

zm =  600*nm
zp = 1000*nm

input = [0., 1.] # [Eix, Eiy]

# Get the S matrix of a left half-infinite block.
Wh, Vh = make_WhVh()


S = zeros(ComplexF64, 4, 4)

S[1:2, 1:2] = aT
S[1:2, 3:4] = aR
S[3:4, 1:2] = bR
S[3:4, 3:4] = bT
