include("./TMM_module.jl")
import .TMM

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

S = zeros(ComplexF64, 4, 4)

# S[1:2, 1:2] = T1
# S[1:2, 3:4] = R1
# S[3:4, 1:2] = R2
# S[3:4, 3:4] = T2

input = [0., 1.] # [Eix, Eiy]
eigvectors, coeff = TMM.left_to_right_field_visualization(dx, dz, Lx, Lz, λ, θ, ϕ, mur, epsr, zm, zp, input)
eigvectors, coeff = TMM.right_to_left_field_visualization(dx, dz, Lx, Lz, λ, θ, ϕ, mur, epsr, zm, zp, input)