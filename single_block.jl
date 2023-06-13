# import Pkg; Pkg.add("NPZ")
# import Pkg; Pkg.add("PythonPlot") # Due to the SSL certificate issue, it cannot be used in LG INNOTEK.
# import Pkg; Pkg.add("PlotlyJS"); Pkg.add("PlotlyBase")
# import Pkg; Pkg.add("UnicodePlots")
using LinearAlgebra
using Printf
using Plots
# gr()
# pythonplot()
# plotlyjs()
# unicodeplots()
using Base
using Distributions
using NPZ

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

function get_HxHy(k0::Number, kzn::Number, A::AbstractMatrix, input::AbstractVector)

    H11 = (1im .* kzn./k0) .* (inv(A) * input)

    Hy = H11[1]
    Hx = H11[2]

    return Hx, Hy
end

function get_Ez(ky_bar, kx_bar, Hy, Hx, epsr)
    Ez = 1im.*(ky_bar.*Hx .- kx_bar.*Hy) ./ epsr
    return Ez
end

function get_Hz(ky_bar, kx_bar, Ey, Ex, mur)
    Hz = 1im.*(ky_bar.*Ex .- kx_bar.*Ey) ./ mur
    return Hz
end

# Eigen vector and eigen value calculation in a single block.
function single_block_eigen(k0, kxn, kyn, mur, epsr)

    kx_bar = kxn / k0
    ky_bar = kyn / k0

    murx = mur
    mury = mur
    murz = mur

    epsrx = epsr
    epsry = epsr
    epsrz = epsr

    A = get_Amatrix(kx_bar, ky_bar, mur, epsr)
    B = get_Bmatrix(kx_bar, ky_bar, mur, epsr)

    C = (k0^2) .* (A * B)

    eigvals, EyEx_eigvecs = eigen(C)

    kzTEp = sqrt(-eigvals[1]) # Left-to-right, TE mode.
    kzTEm =-sqrt(-eigvals[1]) # Right-to-left, TE mode.
    kzTMp = sqrt(-eigvals[2]) # Left-to-right, TM mode.
    kzTMm =-sqrt(-eigvals[2]) # Right-to-left, TM mode.

    EyTE = EyEx_eigvecs[1,1] # Ey, TE mode.
    ExTE = EyEx_eigvecs[2,1] # Ex, TE mode.
    EyTM = EyEx_eigvecs[1,2] # Ey, TM mode.
    ExTM = EyEx_eigvecs[2,2] # Ex, TM mode.

    # Note that H field (magnetic induction vector) is scaled by 1im * sqrt(ε_0/μ_0) for convenience,
    # i.e., H_hat = 1im * sqrt(ε_0/μ_0) * H_bar, and we are calculating H_hat.
    # Also, the order of the eigenvector is reversed, i.e., H11 = [Hy, Hx], not [Hx, Hy].
    HxTEp, HyTEp = get_HxHy(k0, kzTEp, A, [EyTE, ExTE])
    HxTEm, HyTEm = get_HxHy(k0, kzTEm, A, [EyTE, ExTE])
    HxTMp, HyTMp = get_HxHy(k0, kzTMp, A, [EyTM, ExTM])
    HxTMm, HyTMm = get_HxHy(k0, kzTMm, A, [EyTM, ExTM])

    EzTEp = get_Ez(ky_bar, kx_bar, HyTEp, HxTEp, epsrz)
    EzTEm = get_Ez(ky_bar, kx_bar, HyTEm, HxTEm, epsrz)
    EzTMp = get_Ez(ky_bar, kx_bar, HyTMp, HxTMp, epsrz)
    EzTMm = get_Ez(ky_bar, kx_bar, HyTMm, HxTMm, epsrz)

    HzTE = get_Hz(ky_bar, kx_bar, EyTE, ExTE, murz)
    HzTM = get_Hz(ky_bar, kx_bar, EyTM, ExTM, murz)

    eigenvec1 = [ExTE, EyTE, EzTEp, HxTEp, HyTEp, HzTE] # Left-to-right, TE mode.
    eigenvec2 = [ExTE, EyTE, EzTEm, HxTEm, HyTEm, HzTE] # Right-to-left, TE mode.
    eigenvec3 = [ExTM, EyTM, EzTMp, HxTMp, HyTMp, HzTM] # Left-to-right, TM mode.
    eigenvec4 = [ExTM, EyTM, EzTMm, HxTMm, HyTMm, HzTM] # Right-to-left, TM mode.

    _eigvalues = [kzTEp, kzTEm, kzTMp, kzTMm]
    _eigenvecs = hcat(eigenvec1, eigenvec2, eigenvec3, eigenvec4)

    return _eigvalues, _eigenvecs
end

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
function make_WhVh(μ_0, ω::AbstractFloat, kx::AbstractFloat, ky::AbstractFloat, kz::Real)

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

function get_the_left_to_right_operators(μ_0, ω, kx0, ky0, kz0, eigvalues, eigvectors, zm, zp)
    
    Wh, Vh = make_WhVh(μ_0, ω, kx0, ky0, kz0)
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

function get_the_right_to_left_operators(μ_0, ω, kx0, ky0, kz0, eigvalues, eigvectors, zm, zp)
    
    Wh, Vh = make_WhVh(μ_0, ω, kx0, ky0, kz0)
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

"""
    plot_left_to_right_field(input::Tuple{Union{Float64, Int}, Union{Float64, Int}},
    R1::Array{Float64, 2} where size(R1) == (2,2), T1::Array{Float64, 2} where size(R1) == (2,2))

Plot the field profile.

# Arguments
# Returns
# Examples
# Notes
"""
function get_left_to_right_field(dx, dz, Lx, Lz, λ, θ, ϕ, mur, epsr, zm, zp, input::AbstractVector)

    c = 299792458 # m/s.
    k0 = 2*π / λ # wavenumber in free space.
    ω = c*k0
    n = sqrt(mur*epsr)
    μ_0 = 4*π*10^-7
    ε_0 = 8.8541878128e-12
    impedance = sqrt(μ_0 /ε_0)

    # wave vector in free space.
    kx0 = k0 * sin(θ) * cos(ϕ)
    ky0 = k0 * sin(θ) * sin(ϕ)
    kz0 = k0 * cos(θ)

    # wave vector in a single block.
    kxn = k0 * n * sin(θ) * cos(ϕ)
    kyn = k0 * n * sin(θ) * sin(ϕ)

    kzns, eigvectors = single_block_eigen(k0, kxn, kyn, mur, epsr)

    kzTEp = kzns[1]
    kzTEm = kzns[2]
    kzTMp = kzns[3]
    kzTMm = kzns[4]

    eigTEp = eigvectors[:,1]
    eigTEm = eigvectors[:,2]
    eigTMp = eigvectors[:,3]
    eigTMm = eigvectors[:,4]

    @show eigTEp
    @show eigTEm
    @show eigTMp
    @show eigTMm

    # @printf("normalized kx in a single block: %.3f\n", kxn/k0)
    # @printf("normalized ky in a single block: %.3f\n", kyn/k0)
    @printf("normalized kz of TEp in a single block: %.3f\n", kzTEp/k0)
    @printf("normalized kz of TEm in a single block: %.3f\n", kzTEm/k0)
    @printf("normalized kz of TMp in a single block: %.3f\n", kzTMp/k0)
    @printf("normalized kz of TMm in a single block: %.3f\n", kzTMm/k0)

    cap, cam, R1, T1 = get_the_left_to_right_operators(μ_0, ω, kx0, ky0, kz0, kzns, eigvectors, zm, zp)
    cbp, cbm, R2, T2 = get_the_right_to_left_operators(μ_0, ω, kx0, ky0, kz0, kzns, eigvectors, zm, zp)

    # @show cap
    # @show cam
    # @show R1
    # @show T1

    Eix = input[1]
    Eiy = input[2]
    Eiz =-(kx0*Eix + ky0*Eiy) / kz0
    Hix = (ky0*Eiz - kz0*Eiy) / ω / μ_0 # Not H_bar field.
    Hiy = (kz0*Eix - kx0*Eiz) / ω / μ_0
    Hiz = (ky0*Eix - kx0*Eiy) / ω / μ_0

    Erx = dot(R1[2,:], [Eiy, Eix])
    Ery = dot(R1[1,:], [Eiy, Eix])
    Erz = (kx0*Erx + ky0*Ery) / kz0
    Hrx = (ky0*Erz + kz0*Ery) / ω / μ_0
    Hry =-(kz0*Erx + kx0*Erz) / ω / μ_0
    Hrz = (ky0*Erx - kx0*Ery) / ω / μ_0

    Etx = dot(T1[2,:], [Eiy, Eix])
    Ety = dot(T1[1,:], [Eiy, Eix])
    Etz =-(kx0*Etx + ky0*Ety) / kz0
    Htx = (ky0*Etz - kz0*Ety) / ω / μ_0
    Hty = (kz0*Etx - kx0*Etz) / ω / μ_0
    Htz = (ky0*Etx - kx0*Ety) / ω / μ_0

    caTEp = cap[1,:]' * [Eiy, Eix] # The coupling coefficient for left-to-right, TE mode.
    caTMp = cap[2,:]' * [Eiy, Eix] # The coupling coefficient for left-to-right, TM mode.
    caTEm = cam[1,:]' * [Eiy, Eix] # The coupling coefficient for right-to-left, TE mode.
    caTMm = cam[2,:]' * [Eiy, Eix] # The coupling coefficient for right-to-left, TM mode.

    coeff = [caTEp, caTEm, caTMp, caTMm]
    # @show caTEp
    # @show caTMp
    # @show caTEm
    # @show caTMm

    # x = 1000*nm:-10*nm:0 # UnitRange object.
    x = 0:dx:Lx # UnitRange object.
    z = 0:dz:Lz # UnitRange object.
    # x = range(1500, step=-500, stop=0)   # StepRange object.
    # z = range(0, step=500, stop=1500)       # StepRange object.
    y = 0

    # X = repeat(x, 1, length(z))
    # Z = repeat(z, length(x), 1)

    # exit(0)

    leftregion = (z .< zm)
    middleregion = (zm .<= z .< zp)
    rightregion = (zp .<= z)

    # println(x[61])
    # println(middleregion)
    # println(rightregion)

    Ex_tot = zeros(ComplexF64, length(x), length(z))
    Ey_tot = zeros(ComplexF64, length(x), length(z))
    Ez_tot = zeros(ComplexF64, length(x), length(z))
    Hx_tot = zeros(ComplexF64, length(x), length(z))
    Hy_tot = zeros(ComplexF64, length(x), length(z))
    Hz_tot = zeros(ComplexF64, length(x), length(z))

    Ex_inc = zeros(ComplexF64, length(x), length(z))
    Ey_inc = zeros(ComplexF64, length(x), length(z))
    Ez_inc = zeros(ComplexF64, length(x), length(z))
    Hx_inc = zeros(ComplexF64, length(x), length(z))
    Hy_inc = zeros(ComplexF64, length(x), length(z))
    Hz_inc = zeros(ComplexF64, length(x), length(z))

    Ex_ref = zeros(ComplexF64, length(x), length(z))
    Ey_ref = zeros(ComplexF64, length(x), length(z))
    Ez_ref = zeros(ComplexF64, length(x), length(z))
    Hx_ref = zeros(ComplexF64, length(x), length(z))
    Hy_ref = zeros(ComplexF64, length(x), length(z))
    Hz_ref = zeros(ComplexF64, length(x), length(z))

    Ex_trs = zeros(ComplexF64, length(x), length(z))
    Ey_trs = zeros(ComplexF64, length(x), length(z))
    Ez_trs = zeros(ComplexF64, length(x), length(z))
    Hx_trs = zeros(ComplexF64, length(x), length(z))
    Hy_trs = zeros(ComplexF64, length(x), length(z))
    Hz_trs = zeros(ComplexF64, length(x), length(z))

    Ex_TEp = zeros(ComplexF64, length(x), length(z))
    Ey_TEp = zeros(ComplexF64, length(x), length(z))
    Ez_TEp = zeros(ComplexF64, length(x), length(z))
    Hx_TEp = zeros(ComplexF64, length(x), length(z))
    Hy_TEp = zeros(ComplexF64, length(x), length(z))
    Hz_TEp = zeros(ComplexF64, length(x), length(z))

    Ex_TMp = zeros(ComplexF64, length(x), length(z))
    Ey_TMp = zeros(ComplexF64, length(x), length(z))
    Ez_TMp = zeros(ComplexF64, length(x), length(z))
    Hx_TMp = zeros(ComplexF64, length(x), length(z))
    Hy_TMp = zeros(ComplexF64, length(x), length(z))
    Hz_TMp = zeros(ComplexF64, length(x), length(z))

    Ex_TEm = zeros(ComplexF64, length(x), length(z))
    Ey_TEm = zeros(ComplexF64, length(x), length(z))
    Ez_TEm = zeros(ComplexF64, length(x), length(z))
    Hx_TEm = zeros(ComplexF64, length(x), length(z))
    Hy_TEm = zeros(ComplexF64, length(x), length(z))
    Hz_TEm = zeros(ComplexF64, length(x), length(z))

    Ex_TMm = zeros(ComplexF64, length(x), length(z))
    Ey_TMm = zeros(ComplexF64, length(x), length(z))
    Ez_TMm = zeros(ComplexF64, length(x), length(z))
    Hx_TMm = zeros(ComplexF64, length(x), length(z))
    Hy_TMm = zeros(ComplexF64, length(x), length(z))
    Hz_TMm = zeros(ComplexF64, length(x), length(z))

    # For z < z-
    phase_inc = ((kx0.*x) .+ (ky0.*y) .+ (kz0.*(z' .- zm)))[:,leftregion]
    phase_ref = ((kx0.*x) .+ (ky0.*y) .- (kz0.*(z' .- zm)))[:,leftregion]

    Ex_ref[:,leftregion] = (Erx .* exp.(1im .* phase_ref))
    Ey_ref[:,leftregion] = (Ery .* exp.(1im .* phase_ref))
    Ez_ref[:,leftregion] = (Erz .* exp.(1im .* phase_ref))
    Hx_ref[:,leftregion] = (Hrx .* exp.(1im .* phase_ref))
    Hy_ref[:,leftregion] = (Hry .* exp.(1im .* phase_ref))
    Hz_ref[:,leftregion] = (Hrz .* exp.(1im .* phase_ref))

    Ex_inc[:,leftregion] = (Eix .* exp.(1im .* phase_inc))
    Ey_inc[:,leftregion] = (Eiy .* exp.(1im .* phase_inc))
    Ez_inc[:,leftregion] = (Eiz .* exp.(1im .* phase_inc))
    Hx_inc[:,leftregion] = (Hix .* exp.(1im .* phase_inc))
    Hy_inc[:,leftregion] = (Hiy .* exp.(1im .* phase_inc))
    Hz_inc[:,leftregion] = (Hiz .* exp.(1im .* phase_inc))

    # For z- <= z < z+
    phase_TEp = ((kxn.*x) .+ (kyn.*y) .+ (kzTEp.*(z' .- zm)))[:,middleregion]
    phase_TMp = ((kxn.*x) .+ (kyn.*y) .+ (kzTMp.*(z' .- zm)))[:,middleregion]
    phase_TEm = ((kxn.*x) .+ (kyn.*y) .+ (kzTEm.*(z' .- zp)))[:,middleregion]
    phase_TMm = ((kxn.*x) .+ (kyn.*y) .+ (kzTMm.*(z' .- zp)))[:,middleregion]

    Ex_TEp[:,middleregion] = (caTEp .* (eigTEp[1] .* exp.(1im .* phase_TEp)))
    Ex_TEm[:,middleregion] = (caTEm .* (eigTEm[1] .* exp.(1im .* phase_TEm)))
    Ex_TMp[:,middleregion] = (caTMp .* (eigTMp[1] .* exp.(1im .* phase_TMp)))
    Ex_TMm[:,middleregion] = (caTMm .* (eigTMm[1] .* exp.(1im .* phase_TMm)))

    Ey_TEp[:,middleregion] = (caTEp .* (eigTEp[2] .* exp.(1im .* phase_TEp)))
    Ey_TEm[:,middleregion] = (caTEm .* (eigTEm[2] .* exp.(1im .* phase_TEm)))
    Ey_TMp[:,middleregion] = (caTMp .* (eigTMp[2] .* exp.(1im .* phase_TMp)))
    Ey_TMm[:,middleregion] = (caTMm .* (eigTMm[2] .* exp.(1im .* phase_TMm)))

    Ez_TEp[:,middleregion] = (caTEp .* (eigTEp[3] .* exp.(1im .* phase_TEp)))
    Ez_TEm[:,middleregion] = (caTEm .* (eigTEm[3] .* exp.(1im .* phase_TEm)))
    Ez_TMp[:,middleregion] = (caTMp .* (eigTMp[3] .* exp.(1im .* phase_TMp)))
    Ez_TMm[:,middleregion] = (caTMm .* (eigTMm[3] .* exp.(1im .* phase_TMm)))

    Hx_TEp[:,middleregion] = (caTEp .* (eigTEp[4] .* exp.(1im .* phase_TEp)))
    Hx_TEm[:,middleregion] = (caTEm .* (eigTEm[4] .* exp.(1im .* phase_TEm)))
    Hx_TMp[:,middleregion] = (caTMp .* (eigTMp[4] .* exp.(1im .* phase_TMp)))
    Hx_TMm[:,middleregion] = (caTMm .* (eigTMm[4] .* exp.(1im .* phase_TMm)))

    Hy_TEp[:,middleregion] = (caTEp .* (eigTEp[5] .* exp.(1im .* phase_TEp)))
    Hy_TEm[:,middleregion] = (caTEm .* (eigTEm[5] .* exp.(1im .* phase_TEm)))
    Hy_TMp[:,middleregion] = (caTMp .* (eigTMp[5] .* exp.(1im .* phase_TMp)))
    Hy_TMm[:,middleregion] = (caTMm .* (eigTMm[5] .* exp.(1im .* phase_TMm)))

    Hz_TEp[:,middleregion] = (caTEp .* (eigTEp[6] .* exp.(1im .* phase_TEp)))
    Hz_TEm[:,middleregion] = (caTEm .* (eigTEm[6] .* exp.(1im .* phase_TEm)))
    Hz_TMp[:,middleregion] = (caTMp .* (eigTMp[6] .* exp.(1im .* phase_TMp)))
    Hz_TMm[:,middleregion] = (caTMm .* (eigTMm[6] .* exp.(1im .* phase_TMm)))

    # For z+ <= z
    phase_trs = ((kx0.*x) .+ (ky0.*y) .+ (kz0.*(z' .- zp)))[:,rightregion]

    Ex_trs[:,rightregion] = (Etx .* exp.(1im .* phase_trs))
    Ey_trs[:,rightregion] = (Ety .* exp.(1im .* phase_trs))
    Ez_trs[:,rightregion] = (Etz .* exp.(1im .* phase_trs))
    Hx_trs[:,rightregion] = (Htx .* exp.(1im .* phase_trs))
    Hy_trs[:,rightregion] = (Hty .* exp.(1im .* phase_trs))
    Hz_trs[:,rightregion] = (Htz .* exp.(1im .* phase_trs))

    # Add specific fields.
    Ex_inc += Ex_TEp + Ex_TMp + Ex_trs
    Ey_inc += Ey_TEp + Ey_TMp + Ey_trs
    Ez_inc += Ez_TEp + Ez_TMp + Ez_trs
    Hx_inc += Hx_TEp + Hx_TMp + Hx_trs
    Hy_inc += Hy_TEp + Hy_TMp + Hy_trs
    Hz_inc += Hz_TEp + Hz_TMp + Hz_trs

    Ex_ref += Ex_TEm + Ex_TMm
    Ey_ref += Ey_TEm + Ey_TMm
    Ez_ref += Ez_TEm + Ez_TMm
    Hx_ref += Hx_TEm + Hx_TMm
    Hy_ref += Hy_TEm + Hy_TMm
    Hz_ref += Hz_TEm + Hz_TMm

    Ex_tot += Ex_inc + Ex_ref + Ex_TEp + Ex_TEm + Ex_TMp + Ex_TMm + Ex_trs
    Ey_tot += Ey_inc + Ey_ref + Ey_TEp + Ey_TEm + Ey_TMp + Ey_TMm + Ey_trs
    Ez_tot += Ez_inc + Ez_ref + Ez_TEp + Ez_TEm + Ez_TMp + Ez_TMm + Ez_trs
    Hx_tot += Hx_inc + Hx_ref + Hx_TEp + Hx_TEm + Hx_TMp + Hx_TMm + Hx_trs
    Hy_tot += Hy_inc + Hy_ref + Hy_TEp + Hy_TEm + Hy_TMp + Hy_TMm + Hy_trs
    Hz_tot += Hz_inc + Hz_ref + Hz_TEp + Hz_TEm + Hz_TMp + Hz_TMm + Hz_trs

    Hx_inc .*= (1im / impedance)
    Hy_inc .*= (1im / impedance)
    Hz_inc .*= (1im / impedance)
    Hx_ref .*= (1im / impedance)
    Hy_ref .*= (1im / impedance)
    Hz_ref .*= (1im / impedance)
    Hx_tot .*= (1im / impedance)
    Hy_tot .*= (1im / impedance)
    Hz_tot .*= (1im / impedance)

    incfields = cat(real(Ex_inc), real(Hx_inc), real(Ey_inc), real(Hy_inc), real(Ez_inc), real(Hz_inc), dims=3)
    reffields = cat(real(Ex_ref), real(Hx_ref), real(Ey_ref), real(Hy_ref), real(Ez_ref), real(Hz_ref), dims=3)
    trsfields = cat(real(Ex_trs), real(Hx_trs), real(Ey_trs), real(Hy_trs), real(Ez_trs), real(Hz_trs), dims=3)
    totfields = cat(real(Ex_tot), real(Hx_tot), real(Ey_tot), real(Hy_tot), real(Ez_tot), real(Hz_tot), dims=3)

    title=["Ex", "Hx", "Ey", "Hy", "Ez", "Hz"]
    cmins = []
    cmaxs = []

    ps_inc = []
    ps_ref = []
    ps_trs = []
    ps_tot = []
    l = @layout([a d; b e; c f;])

    for (i, slice) in enumerate(eachslice(incfields, dims=3))
        name = title[i]*"_inc"
        cmax = maximum(abs, slice)
        if isapprox(cmax, 0; atol=1e-10)
            p = heatmap(slice, title=name, c=:bwr, clims=(-1,1))
        else
            p = heatmap(slice, title=name, c=:bwr, clims=(-cmax, cmax))
        end
        # gui(p)
        push!(ps_inc, p)
        filename = "./" * name * ".png"
        # println(filename)
        # savefig(p, filename)
    end
    plot(ps_inc..., layout=l, plot_title="inc", size=(1200, 1000))
    savefig("inc.png")
    
    for (i, slice) in enumerate(eachslice(reffields, dims=3))
        name = title[i]*"_ref"
        cmax = maximum(abs, slice)
        if isapprox(cmax, 0; atol=1e-10)
            p = heatmap(slice, title=name, c=:bwr, clims=(-1,1))
        else
            p = heatmap(slice, title=name, c=:bwr, clims=(-cmax, cmax))
        end
        push!(ps_ref, p)
        filename = name * ".png"
        # savefig(filename)
    end
    plot(ps_ref..., layout=l, plot_title="ref", size=(1200, 1000))
    savefig("ref.png")

    for (i, slice) in enumerate(eachslice(trsfields, dims=3))
        name = title[i]*"_trs"
        cmax = maximum(abs, slice)
        if isapprox(cmax, 0; atol=1e-10)
            p = heatmap(slice, title=name, c=:bwr, clims=(-1,1))
        else
            p = heatmap(slice, title=name, c=:bwr, clims=(-cmax, cmax))
        end
        push!(ps_trs, p)
        cmax = maximum(abs, slice)
        push!(cmaxs, cmax)
        filename = name * ".png"
        # savefig(filename)
    end
    plot(ps_trs..., layout=l, plot_title="trs", size=(1200, 1000))
    savefig("trs.png")

    for (i, slice) in enumerate(eachslice(totfields, dims=3))
        name = title[i]*"_tot"
        if isapprox(cmaxs[i], 0; atol=1e-10)
            p = heatmap(slice, title=name, c=:bwr, clims=(-1,1))
        else
            println(cmaxs[i])
            p = heatmap(slice, title=name, c=:bwr, clims=(-cmaxs[i], cmaxs[i]))
        end
        push!(ps_tot, p)
        filename = name * ".png"
        # savefig(filename)
    end
    plot(ps_tot..., layout=l, plot_title="tot", size=(1200, 1000))
    savefig("tot.png")
    # plot(ps..., layout=l, plot_title="trs", size=(1200, 1000), yformatter=:scientific)
    return Ex_tot, Ey_tot, Ez_tot, Hx_tot, Hy_tot, Hz_tot, eigvectors, coeff
end

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
# S[3:4, 1:2] = R1
# S[1:2, 3:4] = R2
# S[3:4, 3:4] = T2

input = [0., 1.] # [Eix, Eiy]
Ex, Ey, Ez, Hx, Hy, Hz, eigvectors, coeff = get_left_to_right_field(dx, dz, Lx, Lz, λ, θ, ϕ, mur, epsr, zm, zp, input)