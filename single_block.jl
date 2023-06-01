using LinearAlgebra
using Printf

um = 1e-6
nm = 1e-9
lambda = 500 * nm
c = 299792458
k0 = 2*π / lambda # wavenumber in free space.
ω = c*k0
mur = 1
epsr = 1.3^2
n = sqrt(mur*epsr)
μ_0 = 4*π*10^-7
ε_0 = 8.8541878128e-12

# polar angle and azimuthal angle.
# θ = π / 12 # radian.
# ϕ = π / 20 # radian.
# θ = π / 6
θ = 0
ϕ = 0
# ϕ = π / 6

# wave vector in free space.
kx0 = k0 * sin(θ) * cos(ϕ)
ky0 = k0 * sin(θ) * sin(ϕ)
kz0 = k0 * cos(θ)

# wave vector in a single block.
kxn = k0 * n * sin(θ) * cos(ϕ)
kyn = k0 * n * sin(θ) * sin(ϕ)

zm = 300*nm
zp = 600*nm

@printf("normalized kx in a single block: %.3f\n", kxn/k0)
@printf("normalized ky in a single block: %.3f\n", kyn/k0)

# Calculate corresponding Ez, Hx, Hy and Hz for the given Ex and Ey.
function get_Amatrix(kx_bar::Number, ky_bar::Number, mur::Number, epsr::Number)

    A = zeros(ComplexF64, 2, 2)

    murx = mur
    mury = mur
    murz = mur

    epsrx = epsr
    epsry = epsr
    epsrz = epsr

    A[1,1] = ky_bar * kx_bar / epsrz
    A[1,2] = murx - ky_bar * ky_bar / epsrz
    A[2,1] =-mury + kx_bar * kx_bar / epsrz
    A[2,2] =-A[1,1]

   return A
end

function get_Bmatrix(kx_bar::Number, ky_bar::Number, mur::Number, epsr::Number)

    B = zeros(ComplexF64, 2, 2)

    murx = mur
    mury = mur
    murz = mur

    epsrx = epsr
    epsry = epsr
    epsrz = epsr

    B[1,1] = ky_bar * kx_bar / murz
    B[1,2] = epsrx - ky_bar * ky_bar / murz
    B[2,1] =-epsry + kx_bar * kx_bar / murz
    B[2,2] =-B[1,1]

   return B
end

# Eigen vector and eigen value calculation in a single block.
function single_block_eigen(k0, kxn, kyn, mur, epsr)

    kx_bar = kxn / k0
    ky_bar = kyn / k0

    A = get_Amatrix(kx_bar, ky_bar, mur, epsr)
    B = get_Bmatrix(kx_bar, ky_bar, mur, epsr)

    C = (k0^2) .* (A * B)

    eigvals, ExEy_eigvecs = eigen(C)

    kz11 = sqrt(-eigvals[1]) # Left-to-right, TE mode.
    kz12 =-sqrt(-eigvals[1]) # Right-to-left, TE mode.
    kz21 = sqrt(-eigvals[2]) # Left-to-right, TM mode.
    kz22 =-sqrt(-eigvals[2]) # Right-to-left, TM mode.

    Ey1 = EyEx_eigvecs[1,1] # Left-to-right, TE mode.
    Ex1 = EyEx_eigvecs[2,1] # Right-to-left, TE mode.
    Ey2 = EyEx_eigvecs[1,2] # Left-to-right, TM mode.
    Ex2 = EyEx_eigvecs[2,2] # Right-to-left, TM mode.

    # Note that H field (magnetic induction vector) is scaled by 1im * sqrt(ε_0/μ_0) for convenience,
    # i.e., H_hat = 1im * sqrt(ε_0/μ_0) * H_bar, and we are calculating H_hat.
    # Also, the order of the eigenvector is reversed, i.e., H11 = [Hy, Hx], not [Hx, Hy].
    H11 = (1im * kz11/k0) .* (inv(A) * [Ey1, Ex1])
    H12 = (1im * kz12/k0) .* (inv(A) * [Ey1, Ex1])
    H21 = (1im * kz21/k0) .* (inv(A) * [Ey2, Ex2])
    H22 = (1im * kz22/k0) .* (inv(A) * [Ey2, Ex2])

    Hx11 = H11[2]
    Hx12 = H12[2]
    Hx21 = H21[2]
    Hx22 = H22[2]

    Hy11 = H11[1]
    Hy12 = H12[1]
    Hy21 = H21[1]
    Hy22 = H22[1]

    Ez11 = 1im*(ky_bar*H11[2] - kx_bar*H11[1]) / epsrz
    Ez12 = 1im*(ky_bar*H12[2] - kx_bar*H12[1]) / epsrz
    Ez21 = 1im*(ky_bar*H21[2] - kx_bar*H21[1]) / epsrz
    Ez22 = 1im*(ky_bar*H22[2] - kx_bar*H22[1]) / epsrz

    Hz1 = 1im*(ky_bar*Ex1 - kx_bar*Ey1) / murz
    Hz2 = 1im*(ky_bar*Ex2 - kx_bar*Ey2) / murz

    eigenvec1 = [Ex1, Ey1, Ez11, Hx11, Hy11, Hz1] # Left-to-right, TE mode.
    eigenvec2 = [Ex1, Ey1, Ez12, Hx12, Hy12, Hz1] # Right-to-left, TE mode.
    eigenvec3 = [Ex2, Ey2, Ez21, Hx21, Hy21, Hz2] # Left-to-right, TM mode.
    eigenvec4 = [Ex2, Ey2, Ez22, Hx22, Hy22, Hz2] # Right-to-left, TM mode.

    # x1_str = join(eigenvec1, ", ")
    # x2_str = join(eigenvec2, ", ")
    # x3_str = join(eigenvec3, ", ")
    # x4_str = join(eigenvec4, ", ")

    # @printf("1st: %.2f, [%s]\n", kz11/k0, x1_str)
    # @printf("2nd: %.2f, [%s]\n", kz12/k0, x2_str)
    # @printf("3rd: %.2f, [%s]\n", kz21/k0, x3_str)
    # @printf("4th: %.2f, [%s]\n", kz22/k0, x4_str)

    # println(typeof(eigenvec1))

    _eigvalues = [kz11, kz12, kz21, kz22]
    _eigenvecs = hcat(eigenvec1, eigenvec2, eigenvec3, eigenvec4)

    return _eigvalues, _eigenvecs
end

kzns, eigvectors = single_block_eigen(k0, kxn, kyn, mur, epsr)

kz1p = kzns[1]
kz1m = kzns[2]
kz2p = kzns[3]
kz2m = kzns[4]

@printf("normalized kz in a single block: %.3f\n", kz1p/k0)
@printf("normalized kz in a single block: %.3f\n", kz1m/k0)
@printf("normalized kz in a single block: %.3f\n", kz2p/k0)
@printf("normalized kz in a single block: %.3f\n", kz2m/k0)

"""
    make_WhVh(ω, kx, ky, kz, eigvecs, z)

Make Wh and Vh matrices.

# Arguments
- 'ω::Float': angular frequency of the input wave.
- 'kx::Float': kx in free space.
- 'ky::Float': ky in free space.
- 'kz::Float': kz in free space.

# Returns

# Examples
'''julia
julia>
'''

# Notes
"""
function make_WhVh(ω::AbstractFloat, kx::AbstractFloat, ky::AbstractFloat, kz::Real)

    Wh = I # The most Julianic way of expressing the identity matrix.
    Vh = zeros(ComplexF64, 2, 2)

    Vh[1,1] = kx * ky / kz / ω / μ_0
    Vh[1,2] = (kx^2 + kz^2) / kz / ω / μ_0
    Vh[2,1] =-(ky^2 + kz^2) / kz / ω / μ_0
    Vh[2,2] =-kx * ky / kz / ω / μ_0

    return Wh, Vh
end

"""
    make_WpVp(eigvecs, z)

Make some matrices and calculate the coupling coeffient matrix operators.

# Arguments

# Returns

# Examples
```julia
julia>
```

# Notes
"""
function make_WpVp(eigvalues::AbstractVector, eigvecs::AbstractMatrix, z::Real)

    kz1p = eigvalues[1]
    kz1m = eigvalues[2]
    kz2p = eigvalues[3]
    kz2m = eigvalues[4]

    ExTEp = eigvecs[1,1]
    ExTEm = eigvecs[1,2]
    ExTMp = eigvecs[1,3]
    ExTMm = eigvecs[1,4]

    EyTEp = eigvecs[2,1]
    EyTEm = eigvecs[2,2]
    EyTMp = eigvecs[2,3]
    EyTMm = eigvecs[2,4]

    EzTEp = eigvecs[3,1]
    EzTEm = eigvecs[3,2]
    EzTMp = eigvecs[3,3]
    EzTMm = eigvecs[3,4]

    HxTEp = eigvecs[4,1]
    HxTEm = eigvecs[4,2]
    HxTMp = eigvecs[4,3]
    HxTMm = eigvecs[4,4]

    HyTEp = eigvecs[5,1]
    HyTEm = eigvecs[5,2]
    HyTMp = eigvecs[5,3]
    HyTMm = eigvecs[5,4]

    HzTEp = eigvecs[6,1]
    HzTEm = eigvecs[6,2]
    HzTMp = eigvecs[6,3]
    HzTMm = eigvecs[6,4]

    Wp = zeros(ComplexF64, 2, 2)
    Vp = zeros(ComplexF64, 2, 2)

    Wp[1,1] = EyTEp * exp(1im*kz1p*z)
    Wp[1,2] = EyTMp * exp(1im*kz2p*z)
    Wp[2,1] = ExTEp * exp(1im*kz1p*z)
    Wp[2,2] = ExTMp * exp(1im*kz2p*z)

    Vp[1,1] = HyTEp * exp(1im*kz1p*z)
    Vp[1,2] = HyTMp * exp(1im*kz2p*z)
    Vp[2,1] = HxTEp * exp(1im*kz1p*z)
    Vp[2,2] = HxTMp * exp(1im*kz2p*z)

    return Wp, Vp
end

function make_WmVm(eigvalues::AbstractVector, eigvecs::AbstractMatrix, z::Real)

    kz1p = eigvalues[1]
    kz1m = eigvalues[2]
    kz2p = eigvalues[3]
    kz2m = eigvalues[4]

    ExTEp = eigvecs[1,1]
    ExTEm = eigvecs[1,2]
    ExTMp = eigvecs[1,3]
    ExTMm = eigvecs[1,4]

    EyTEp = eigvecs[2,1]
    EyTEm = eigvecs[2,2]
    EyTMp = eigvecs[2,3]
    EyTMm = eigvecs[2,4]

    EzTEp = eigvecs[3,1]
    EzTEm = eigvecs[3,2]
    EzTMp = eigvecs[3,3]
    EzTMm = eigvecs[3,4]

    HxTEp = eigvecs[4,1]
    HxTEm = eigvecs[4,2]
    HxTMp = eigvecs[4,3]
    HxTMm = eigvecs[4,4]

    HyTEp = eigvecs[5,1]
    HyTEm = eigvecs[5,2]
    HyTMp = eigvecs[5,3]
    HyTMm = eigvecs[5,4]

    HzTEp = eigvecs[6,1]
    HzTEm = eigvecs[6,2]
    HzTMp = eigvecs[6,3]
    HzTMm = eigvecs[6,4]

    Wm = zeros(ComplexF64, 2, 2)
    Vm = zeros(ComplexF64, 2, 2)

    Wm[1,1] = EyTEm * exp(1im*kz1m*z)
    Wm[1,2] = EyTMm * exp(1im*kz2m*z)
    Wm[2,1] = ExTEm * exp(1im*kz1m*z)
    Wm[2,2] = ExTMm * exp(1im*kz2m*z)

    Vm[1,1] = HyTEm * exp(1im*kz1m*z)
    Vm[1,2] = HyTMm * exp(1im*kz2m*z)
    Vm[2,1] = HxTEm * exp(1im*kz1m*z)
    Vm[2,2] = HxTMm * exp(1im*kz2m*z)

    return Wm, Vm
end

function get_the_left_to_right_operators(ω, kx0, ky0, kz0, eigvalues, eigvectors, zm, zp)
    
    Wh, Vh = make_WhVh(ω, kx0, ky0, kz0)
    Wp0, Vp0 = make_WpVp(eigvalues, eigvectors, 0.)
    Wpp, Vpp = make_WpVp(eigvalues, eigvectors, zp-zm)
    Wm0, Vm0 = make_WmVm(eigvalues, eigvectors, 0)
    Wmm, Vmm = make_WmVm(eigvalues, eigvectors, zm-zp)

    A = inv(Wh)*Wp0 + inv(Vh)*Vp0
    B = inv(Wh)*Wmm + inv(Vh)*Vmm
    C = inv(Wh)*Wpp - inv(Vh)*Vpp
    D = inv(Wh)*Wm0 - inv(Vh)*Vm0

    E = zeros(ComplexF64, 4, 4)
    E[1:2,1:2] = A
    E[1:2,3:4] = B
    E[3:4,1:2] = C
    E[3:4,3:4] = D

    u = zeros(Float64, 4, 2)
    u[1,1] = 2
    u[2,2] = 2

    ca = inv(E) * u

    cap = ca[1:2,1:2]
    cam = ca[3:4,1:2]

    R = inv(Wh) * (Wp0*cap + Wmm*cam - Wh*I)
    T = inv(Wh) * (Wpp*cap + Wm0*cam)

    return cap, cam, R, T
end

function get_the_right_to_left_operators(ω, kx0, ky0, kz0, eigvalues, eigvectors, zm, zp)
    
    Wh, Vh = make_WhVh(ω, kx0, ky0, kz0)
    Wp0, Vp0 = make_WpVp(eigvalues, eigvectors, 0.)
    Wpp, Vpp = make_WpVp(eigvalues, eigvectors, zp-zm)
    Wm0, Vm0 = make_WmVm(eigvalues, eigvectors, 0)
    Wmm, Vmm = make_WmVm(eigvalues, eigvectors, zm-zp)

    A = inv(Wh)*Wp0 + inv(Vh)*Vp0
    B = inv(Wh)*Wmm + inv(Vh)*Vmm
    C = inv(Wh)*Wpp - inv(Vh)*Vpp
    D = inv(Wh)*Wm0 - inv(Vh)*Vm0

    E = zeros(ComplexF64, 4, 4)
    E[1:2,1:2] = A
    E[1:2,3:4] = B
    E[3:4,1:2] = C
    E[3:4,3:4] = D

    u = zeros(Float64, 4, 2)
    u[3,1] = 2
    u[4,2] = 2

    cb = inv(E) * u

    cbp = cb[1:2,1:2]
    cbm = cb[3:4,1:2]

    R = inv(Wh) * (Wp0*cbp + Wm0*cbm - Wh*I)
    T = inv(Wh) * (Wp0*cbp + Wmm*cbm)

    return cbp, cbm, R, T
end

cap, cam, R1, T1 = get_the_left_to_right_operators(ω, kx0, ky0, kz0, kzns, eigvectors, zm, zp)
cbp, cbm, R2, T2 = get_the_right_to_left_operators(ω, kx0, ky0, kz0, kzns, eigvectors, zm, zp)

S = zeros(ComplexF64, 4, 4)

S[1:2, 1:2] = T1
S[3:4, 1:2] = R1
S[1:2, 3:4] = R2
S[3:4, 3:4] = T2

"""
    plot_left_to_right_field(input::Tuple{Union{Float64, Int}, Union{Float64, Int}},
    R1::Array{Float64, 2} where size(R1) == (2,2), T1::Array{Float64, 2} where size(R1) == (2,2))

Plot the field profile.

# Arguments
# Returns
# Examples
# Notes
"""
function plot_left_to_right_field(ω, μ_0, kx0, ky0, kz0, zm, zp, kzns,
    input::AbstractVector, eigvectors::AbstractVector,
    cap::AbstractMatrix, cam::AbstractMatrix,
    R1::AbstractMatrix{Float64}, T1::AbstractMatrix)

    kzTEp = kzns[1]
    kzTEm = kzns[2]
    kzTMp = kzns[3]
    kzTMm = kzns[4]

    eigTEp = eigvectors[:,1]
    eigTEm = eigvectors[:,2]
    eigTMp = eigvectors[:,3]
    eigTMm = eigvectors[:,4]

    Eiy = input[1]
    Eix = input[2]
    Eiz =-(kx0*Eix + ky0*Eiy) / kz0
    Hiy = (kz0*Eix - kx0*Eiz) / ω / μ_0
    Hix = (ky0*Eiz - kz0*Eiy) / ω / μ_0
    Hiz = (ky0*Eix - kx0*Eiy) / ω / μ_0

    Ery = R1[1,:] * input
    Erx = R1[2,:] * input
    Erz = (kx0*Erx + ky0*Ery) / kz0
    Hry =-(kz0*Erx + kx0*Erz) / ω / μ_0
    Hrx = (ky0*Erz + kz0*Ery) / ω / μ_0
    Hrz = (ky0*Erx - kx0*Ery) / ω / μ_0

    Ety = T1[1,:] * input
    Etx = T1[2,:] * input
    Etz =-(kx0*Etx + ky0*Ety) / kz0
    Hty = (kz0*Etx - kx0*Etz) / ω / μ_0
    Htx = (ky0*Etz - kz0*Ety) / ω / μ_0
    Htz = (ky0*Etx - kx0*Ety) / ω / μ_0

    capp = cap * input # The coupling coefficient for left-to-right, {TE mode, TM mode}.
    camm = cam * input # The coupling coefficient for right-to-left, {TE mode, TM mode}.

    x = 0:1*nm:1500*nm
    z = 0:1*nm:1500*nm
    y = 0

    Z, X = ndgrid(z, x, indexing=:xy)
    field = zero(Z)

    # For z < z-
    phase_inc = (kx0 .* X) .+ (ky0 * y) + (kz0 .* (Z[Z .< zm] .- zm))
    phase_ref = (kx0 .* X) .+ (ky0 * y) - (kz0 .* (Z[Z .< zm] .- zm))
    leftfield = ([Eix, Eiy, Eiz, Hix, Hiy, Hiz] .* exp(1im .* phase_inc)) + 
                ([Erx, Ery, Erz, Hrx, Hry, Hrz] .* exp(1im .* phase_ref))

    # For z- <= z < z+
    phase_TEp = (kxn .* X) .+ (kyn * y) + (kzTEp .* (Z[(Z .=> zm) & (Z .< zp)] .- zm))
    phase_TMp = (kxn .* X) .+ (kyn * y) + (kzTMp .* (Z .- zm))
    phase_TEm = (kxn .* X) .+ (kyn * y) + (kzTEm .* (Z .- zp))
    phase_TMm = (kxn .* X) .+ (kyn * y) + (kzTMm .* (Z .- zp))
    middlefield = (capp[1] .* (eigTEp .* phase_TEp)) + (camm[1] .* (eigTEm .* phase_TEm)) +
                  (capp[2] .* (eigTMp .* phase_TMp)) + (camm[2] .* (eigTMm .* phase_TMm))

    # For z+ <= z
    phase_trs = (kx0 .* X) .+ (ky0 * y) + (kz0 .* (Z .- zp))
    rightfield = ([Etx, Ety, Etz, Htx, Hty, Htz] .* exp(1im .* phase_trs))

    return
end