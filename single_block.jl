using Pkg

#Pkg.add("Plots")
#Pkg.add("NPZ")
#Pkg.add("LinearAlgebra")
#Pkg.add("Distributions")
#Pkg.add("PythonPlot") # Due to the SSL certificate issue, it cannot be used in LG INNOTEK.
#Pkg.add("PlotlyJS"); Pkg.add("PlotlyBase")
#Pkg.add("UnicodePlots")

"""
# Check if the package is installed
if "Plots" in keys(Pkg.dependencies())
    println("The Plots package is installed.")
else
    println("The Plots package is not installed.")
	Pkg.add("Plots")
end

if "NPZ" in keys(Pkg.dependencies())
    println("The NPZ package is installed.")
else
    println("The NPZ package is not installed.")
	Pkg.add("NPZ")
end

if "LinearAlgebra" in keys(Pkg.dependencies())
    println("The LinearAlgebra package is installed.")
else
    println("The LinearAlgebra package is not installed.")
	Pkg.add("LinearAlgebra")
end

if "Distributions" in keys(Pkg.dependencies())
    println("The Distributions package is installed.")
else
    println("The Distributions package is not installed.")
	Pkg.add("Distributions")
end
"""

include("./TMM_module.jl")

using Base
using Distributions
using LinearAlgebra
using NPZ
using Printf
using Plots
using .single_block

# Choose the backends of the Plots package.
gr()
# pythonplot()
# plotlyjs()
# unicodeplots()


"""
    plot_left_to_right_field(input::Tuple{Union{Float64, Int}, Union{Float64, Int}},
    R1::Array{Float64, 2} where size(R1) == (2,2), T1::Array{Float64, 2} where size(R1) == (2,2))

Plot the field profile.

# Arguments
# Returns
# Examples
# Notes
"""
function left_to_right_field_visualization(dx, dz, Lx, Lz, λ, θ, ϕ, mur, epsr, zm, zp, input::AbstractVector)

    c = 299792458 # m/s.
    k0 = 2*π / λ # wavenumber in free space.
    ω = c*k0
    n = sqrt(mur*epsr)
    μ_0 = 4*π*10^-7
    ε_0 = 8.8541878128e-12
    impedance = sqrt(μ_0 /ε_0)

    # wave vector in free space.
    # Note the the tangential component of the wavevector at the boundary
    # is continuous, which means that they only differ in the refractive index.
    kx0 = k0 * sin(θ) * cos(ϕ)
    ky0 = k0 * sin(θ) * sin(ϕ)
    kz0 = k0 * cos(θ)

    # x = 1000*nm:-10*nm:0 # UnitRange object.
    x = 0:dx:Lx # UnitRange object.
    z = 0:dz:Lz # UnitRange object.
    # x = range(1500, step=-500, stop=0)   # StepRange object.
    # z = range(0, step=500, stop=1500)       # StepRange object.
    y = 0
    # X = repeat(x, 1, length(z))
    # Z = repeat(z, length(x), 1)
    # exit(0)

    kzns, eigvectors = get_eigenvectors(k0, kx0, ky0, mur, epsr, impedance)

    kzTEp = kzns[1]
    kzTEm = kzns[2]
    kzTMp = kzns[3]
    kzTMm = kzns[4]

    eigTEp = eigvectors[:,1]
    eigTEm = eigvectors[:,2]
    eigTMp = eigvectors[:,3]
    eigTMm = eigvectors[:,4]

    # @show eigTEp
    # @show eigTEm
    # @show eigTMp
    # @show eigTMm

    kx_bar = kx0 / k0
    ky_bar = ky0 / k0
    kzTEp_bar = kzTEp / k0
    kzTEm_bar = kzTEm / k0
    kzTMp_bar = kzTMp / k0
    kzTMm_bar = kzTMm / k0
    norm_mag = sqrt(kx_bar^2 + ky_bar^2 + kzTEp_bar^2)

    @printf("normalized kx in a single block: %.3f\n", kx_bar)
    @printf("normalized ky in a single block: %.3f\n", ky_bar)
    @printf("normalized kz of TEp in a single block: %.3f\n", kzTEp_bar)
    @printf("normalized kz of TEm in a single block: %.3f\n", kzTEm_bar)
    @printf("normalized kz of TMp in a single block: %.3f\n", kzTMp_bar)
    @printf("normalized kz of TMm in a single block: %.3f\n", kzTMm_bar)
    @printf("normalized magnitude of the wavevector in a single block: %.3f\n", norm_mag)

    cap, cam, R1, T1 = get_the_left_to_right_operators(μ_0, ω, kx0, ky0, kz0, kzns, eigvectors, zm, zp)

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
    phase_TEp = ((kx0.*x) .+ (ky0.*y) .+ (kzTEp.*(z' .- zm)))[:,middleregion]
    phase_TMp = ((kx0.*x) .+ (ky0.*y) .+ (kzTMp.*(z' .- zm)))[:,middleregion]
    phase_TEm = ((kx0.*x) .+ (ky0.*y) .+ (kzTEm.*(z' .- zp)))[:,middleregion]
    phase_TMm = ((kx0.*x) .+ (ky0.*y) .+ (kzTMm.*(z' .- zp)))[:,middleregion]

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
    savefig("LtoR_inc.png")
    
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
    savefig("LtoR_ref.png")

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
    savefig("LtoR_trs.png")

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
    savefig("LtoR_tot.png")
    # plot(ps..., layout=l, plot_title="trs", size=(1200, 1000), yformatter=:scientific)
    
   return eigvectors, coeff
end

function right_to_left_field_visualization(dx, dz, Lx, Lz, λ, θ, ϕ, mur, epsr, zm, zp, input::AbstractVector)

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

    # x = 1000*nm:-10*nm:0 # UnitRange object.
    x = 0:dx:Lx # UnitRange object.
    z = 0:dz:Lz # UnitRange object.
    # x = range(1500, step=-500, stop=0)   # StepRange object.
    # z = range(0, step=500, stop=1500)       # StepRange object.
    y = 0
    # X = repeat(x, 1, length(z))
    # Z = repeat(z, length(x), 1)
    # exit(0)

    kzns, eigvectors = get_eigenvectors(k0, kx0, ky0, mur, epsr, impedance)

    kzTEp = kzns[1]
    kzTEm = kzns[2]
    kzTMp = kzns[3]
    kzTMm = kzns[4]

    eigTEp = eigvectors[:,1]
    eigTEm = eigvectors[:,2]
    eigTMp = eigvectors[:,3]
    eigTMm = eigvectors[:,4]

    cbp, cbm, R, T = get_the_right_to_left_operators(μ_0, ω, kx0, ky0, kz0, kzns, eigvectors, zm, zp)

    Eix = input[1]
    Eiy = input[2]
    Eiz = (kx0*Eix + ky0*Eiy) / kz0
    Hix = (ky0*Eiz + kz0*Eiy) / ω / μ_0 # Not H_bar field.
    Hiy =-(kz0*Eix + kx0*Eiz) / ω / μ_0
    Hiz = (ky0*Eix + kx0*Eiy) / ω / μ_0

    Erx = dot(R[2,:], [Eiy, Eix])
    Ery = dot(R[1,:], [Eiy, Eix])
    Erz =-(kx0*Erx + ky0*Ery) / kz0
    Hrx = (ky0*Erz - kz0*Ery) / ω / μ_0
    Hry = (kz0*Erx - kx0*Erz) / ω / μ_0
    Hrz = (ky0*Erx + kx0*Ery) / ω / μ_0

    Etx = dot(T[2,:], [Eiy, Eix])
    Ety = dot(T[1,:], [Eiy, Eix])
    Etz =-(kx0*Etx + ky0*Ety) / kz0
    Htx = (ky0*Etz + kz0*Ety) / ω / μ_0
    Hty =-(kz0*Etx + kx0*Etz) / ω / μ_0
    Htz = (ky0*Etx + kx0*Ety) / ω / μ_0

    cbTEp = cbp[1,:]' * [Eiy, Eix] # The coupling coefficient for left-to-right, TE mode.
    cbTMp = cbp[2,:]' * [Eiy, Eix] # The coupling coefficient for left-to-right, TM mode.
    cbTEm = cbm[1,:]' * [Eiy, Eix] # The coupling coefficient for right-to-left, TE mode.
    cbTMm = cbm[2,:]' * [Eiy, Eix] # The coupling coefficient for right-to-left, TM mode.

    coeff = [cbTEp, cbTEm, cbTMp, cbTMm]

    leftregion = (z .< zm)
    middleregion = (zm .<= z .< zp)
    rightregion = (zp .<= z)

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
    phase_trs = ((kx0.*x) .+ (ky0.*y) .- (kz0.*(z' .- zm)))[:,leftregion]

    Ex_trs[:,leftregion] = (Etx .* exp.(1im .* phase_trs))
    Ey_trs[:,leftregion] = (Ety .* exp.(1im .* phase_trs))
    Ez_trs[:,leftregion] = (Etz .* exp.(1im .* phase_trs))
    Hx_trs[:,leftregion] = (Htx .* exp.(1im .* phase_trs))
    Hy_trs[:,leftregion] = (Hty .* exp.(1im .* phase_trs))
    Hz_trs[:,leftregion] = (Htz .* exp.(1im .* phase_trs))

    # For z- <= z < z+
    phase_TEp = ((kx0.*x) .+ (ky0.*y) .+ (kzTEp.*(z' .- zm)))[:,middleregion]
    phase_TMp = ((kx0.*x) .+ (ky0.*y) .+ (kzTMp.*(z' .- zm)))[:,middleregion]
    phase_TEm = ((kx0.*x) .+ (ky0.*y) .+ (kzTEm.*(z' .- zp)))[:,middleregion]
    phase_TMm = ((kx0.*x) .+ (ky0.*y) .+ (kzTMm.*(z' .- zp)))[:,middleregion]

    Ex_TEp[:,middleregion] = (cbTEp .* (eigTEp[1] .* exp.(1im .* phase_TEp)))
    Ex_TEm[:,middleregion] = (cbTEm .* (eigTEm[1] .* exp.(1im .* phase_TEm)))
    Ex_TMp[:,middleregion] = (cbTMp .* (eigTMp[1] .* exp.(1im .* phase_TMp)))
    Ex_TMm[:,middleregion] = (cbTMm .* (eigTMm[1] .* exp.(1im .* phase_TMm)))

    Ey_TEp[:,middleregion] = (cbTEp .* (eigTEp[2] .* exp.(1im .* phase_TEp)))
    Ey_TEm[:,middleregion] = (cbTEm .* (eigTEm[2] .* exp.(1im .* phase_TEm)))
    Ey_TMp[:,middleregion] = (cbTMp .* (eigTMp[2] .* exp.(1im .* phase_TMp)))
    Ey_TMm[:,middleregion] = (cbTMm .* (eigTMm[2] .* exp.(1im .* phase_TMm)))

    Ez_TEp[:,middleregion] = (cbTEp .* (eigTEp[3] .* exp.(1im .* phase_TEp)))
    Ez_TEm[:,middleregion] = (cbTEm .* (eigTEm[3] .* exp.(1im .* phase_TEm)))
    Ez_TMp[:,middleregion] = (cbTMp .* (eigTMp[3] .* exp.(1im .* phase_TMp)))
    Ez_TMm[:,middleregion] = (cbTMm .* (eigTMm[3] .* exp.(1im .* phase_TMm)))

    Hx_TEp[:,middleregion] = (cbTEp .* (eigTEp[4] .* exp.(1im .* phase_TEp)))
    Hx_TEm[:,middleregion] = (cbTEm .* (eigTEm[4] .* exp.(1im .* phase_TEm)))
    Hx_TMp[:,middleregion] = (cbTMp .* (eigTMp[4] .* exp.(1im .* phase_TMp)))
    Hx_TMm[:,middleregion] = (cbTMm .* (eigTMm[4] .* exp.(1im .* phase_TMm)))

    Hy_TEp[:,middleregion] = (cbTEp .* (eigTEp[5] .* exp.(1im .* phase_TEp)))
    Hy_TEm[:,middleregion] = (cbTEm .* (eigTEm[5] .* exp.(1im .* phase_TEm)))
    Hy_TMp[:,middleregion] = (cbTMp .* (eigTMp[5] .* exp.(1im .* phase_TMp)))
    Hy_TMm[:,middleregion] = (cbTMm .* (eigTMm[5] .* exp.(1im .* phase_TMm)))

    Hz_TEp[:,middleregion] = (cbTEp .* (eigTEp[6] .* exp.(1im .* phase_TEp)))
    Hz_TEm[:,middleregion] = (cbTEm .* (eigTEm[6] .* exp.(1im .* phase_TEm)))
    Hz_TMp[:,middleregion] = (cbTMp .* (eigTMp[6] .* exp.(1im .* phase_TMp)))
    Hz_TMm[:,middleregion] = (cbTMm .* (eigTMm[6] .* exp.(1im .* phase_TMm)))

    # For z+ <= z
    phase_inc = ((kx0.*x) .+ (ky0.*y) .- (kz0.*(z' .- zp)))[:,rightregion]
    phase_ref = ((kx0.*x) .+ (ky0.*y) .+ (kz0.*(z' .- zp)))[:,rightregion]

    Ex_inc[:,rightregion] = (Eix .* exp.(1im .* phase_inc))
    Ey_inc[:,rightregion] = (Eiy .* exp.(1im .* phase_inc))
    Ez_inc[:,rightregion] = (Eiz .* exp.(1im .* phase_inc))
    Hx_inc[:,rightregion] = (Hix .* exp.(1im .* phase_inc))
    Hy_inc[:,rightregion] = (Hiy .* exp.(1im .* phase_inc))
    Hz_inc[:,rightregion] = (Hiz .* exp.(1im .* phase_inc))

    Ex_ref[:,rightregion] = (Erx .* exp.(1im .* phase_ref))
    Ey_ref[:,rightregion] = (Ery .* exp.(1im .* phase_ref))
    Ez_ref[:,rightregion] = (Erz .* exp.(1im .* phase_ref))
    Hx_ref[:,rightregion] = (Hrx .* exp.(1im .* phase_ref))
    Hy_ref[:,rightregion] = (Hry .* exp.(1im .* phase_ref))
    Hz_ref[:,rightregion] = (Hrz .* exp.(1im .* phase_ref))

    # Add specific fields.
    Ex_inc += Ex_TEm + Ex_TMm + Ex_trs
    Ey_inc += Ey_TEm + Ey_TMm + Ey_trs
    Ez_inc += Ez_TEm + Ez_TMm + Ez_trs
    Hx_inc += Hx_TEm + Hx_TMm + Hx_trs
    Hy_inc += Hy_TEm + Hy_TMm + Hy_trs
    Hz_inc += Hz_TEm + Hz_TMm + Hz_trs

    Ex_ref += Ex_TEp + Ex_TMp
    Ey_ref += Ey_TEp + Ey_TMp
    Ez_ref += Ez_TEp + Ez_TMp
    Hx_ref += Hx_TEp + Hx_TMp
    Hy_ref += Hy_TEp + Hy_TMp
    Hz_ref += Hz_TEp + Hz_TMp

    Ex_tot += Ex_inc + Ex_ref + Ex_TEp + Ex_TEm + Ex_TMp + Ex_TMm + Ex_trs
    Ey_tot += Ey_inc + Ey_ref + Ey_TEp + Ey_TEm + Ey_TMp + Ey_TMm + Ey_trs
    Ez_tot += Ez_inc + Ez_ref + Ez_TEp + Ez_TEm + Ez_TMp + Ez_TMm + Ez_trs
    Hx_tot += Hx_inc + Hx_ref + Hx_TEp + Hx_TEm + Hx_TMp + Hx_TMm + Hx_trs
    Hy_tot += Hy_inc + Hy_ref + Hy_TEp + Hy_TEm + Hy_TMp + Hy_TMm + Hy_trs
    Hz_tot += Hz_inc + Hz_ref + Hz_TEp + Hz_TEm + Hz_TMp + Hz_TMm + Hz_trs

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
    end
    plot(ps_inc..., layout=l, plot_title="inc", size=(1200, 1000))
    savefig("RtoL_inc.png")
    
    for (i, slice) in enumerate(eachslice(reffields, dims=3))
        name = title[i]*"_ref"
        cmax = maximum(abs, slice)
        if isapprox(cmax, 0; atol=1e-10)
            p = heatmap(slice, title=name, c=:bwr, clims=(-1,1))
        else
            p = heatmap(slice, title=name, c=:bwr, clims=(-cmax, cmax))
        end
        push!(ps_ref, p)
    end
    plot(ps_ref..., layout=l, plot_title="ref", size=(1200, 1000))
    savefig("RtoL_ref.png")

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
    end
    plot(ps_trs..., layout=l, plot_title="trs", size=(1200, 1000))
    savefig("RtoL_trs.png")

    for (i, slice) in enumerate(eachslice(totfields, dims=3))
        name = title[i]*"_tot"
        if isapprox(cmaxs[i], 0; atol=1e-10)
            p = heatmap(slice, title=name, c=:bwr, clims=(-1,1))
        else
            println(cmaxs[i])
            p = heatmap(slice, title=name, c=:bwr, clims=(-cmaxs[i], cmaxs[i]))
        end
        push!(ps_tot, p)
    end
    plot(ps_tot..., layout=l, plot_title="tot", size=(1200, 1000))
    savefig("RtoL_tot.png")
    # plot(ps..., layout=l, plot_title="trs", size=(1200, 1000), yformatter=:scientific)
    
   return eigvectors, coeff
end

using .single_block


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
eigvectors, coeff = left_to_right_field_visualization(dx, dz, Lx, Lz, λ, θ, ϕ, mur, epsr, zm, zp, input)
eigvectors, coeff = right_to_left_field_visualization(dx, dz, Lx, Lz, λ, θ, ϕ, mur, epsr, zm, zp, input)