module TMM

	# The following packages are required in this module.
	# using LinearAlgebra
	import LinearAlgebra as la

	export _get_Amatrix, _get_Bmatrix
	export _get_HxHy, _get_Ez, _get_Hz 
	export _make_WhVh, _make_WpVp, _make_WmVm
	export get_eigenvectors
	export get_the_left_to_right_operators
	export get_the_right_to_left_operators
	export left_to_right_field_visualization
	export right_to_left_field_visualization
	export get_scc_of_a_block
	export get_scc_of_each_n_blocks
	export redheffer
	export redheffer_n_blocks
	export redheffer_two_blocks
	export redheffer_combine_all_blocks

    # Calculate corresponding Ez, Hx, Hy and Hz for the given Ex and Ey.
	function _get_Amatrix(
		kx_bar::Number, 
		ky_bar::Number, 
		mur::Number, 
		epsr::Number
	)

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

	function _get_Bmatrix(
		kx_bar::Number, 
		ky_bar::Number, 
		mur::Number, 
		epsr::Number
	)

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

	function _get_HxHy(
		k0::Number, 
		kzn::Number, 
		A::AbstractMatrix, 
		input::AbstractVector, 
		impedance::Number)

		H11 = (1im .* kzn./k0) .* (inv(A) * input)

		Hy = H11[1]
		Hx = H11[2]

		Hx *= (1im / impedance)
		Hy *= (1im / impedance)

		return Hx, Hy
	end

	function _get_Ez(ky_bar, kx_bar, Hy, Hx, epsr, impedence)
		Hy *= (impedence / 1im)
		Hx *= (impedence / 1im)
		Ez = 1im.*(ky_bar.*Hx .- kx_bar.*Hy) ./ epsr
		return Ez
	end

	function _get_Hz(ky_bar, kx_bar, Ey, Ex, mur, impedance)
		Hz = 1im.*(ky_bar.*Ex .- kx_bar.*Ey) ./ mur
		Hz*= (1im / impedance)
		return Hz
	end

    # Obtain the eigen vectors and eigenvalues in a medium.
	function get_eigenvectors(k0, kx0, ky0, mur, epsr)

		μ_0 = 4*π*10^-7
		ε_0 = 8.8541878128e-12
		impedance = sqrt(μ_0 /ε_0)

		kx_bar = kx0 / k0 # it is equal to sin(θ) cos(ϕ)
		ky_bar = ky0 / k0 # it is equal to sin(θ) sin(ϕ)

		murx = mur
		mury = mur
		murz = mur

		epsrx = epsr
		epsry = epsr
		epsrz = epsr

		A = _get_Amatrix(kx_bar, ky_bar, mur, epsr)
		B = _get_Bmatrix(kx_bar, ky_bar, mur, epsr)

		C = (k0^2) .* (A * B)

		eigvals, EyEx_eigvecs = la.eigen(C)

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
		HxTEp, HyTEp = _get_HxHy(k0, kzTEp, A, [EyTE, ExTE], impedance)
		HxTEm, HyTEm = _get_HxHy(k0, kzTEm, A, [EyTE, ExTE], impedance)
		HxTMp, HyTMp = _get_HxHy(k0, kzTMp, A, [EyTM, ExTM], impedance)
		HxTMm, HyTMm = _get_HxHy(k0, kzTMm, A, [EyTM, ExTM], impedance)

		EzTEp = _get_Ez(ky_bar, kx_bar, HyTEp, HxTEp, epsrz, impedance)
		EzTEm = _get_Ez(ky_bar, kx_bar, HyTEm, HxTEm, epsrz, impedance)
		EzTMp = _get_Ez(ky_bar, kx_bar, HyTMp, HxTMp, epsrz, impedance)
		EzTMm = _get_Ez(ky_bar, kx_bar, HyTMm, HxTMm, epsrz, impedance)

		HzTE = _get_Hz(ky_bar, kx_bar, EyTE, ExTE, murz, impedance)
		HzTM = _get_Hz(ky_bar, kx_bar, EyTM, ExTM, murz, impedance)

		eigenvec1 = [ExTE, EyTE, EzTEp, HxTEp, HyTEp, HzTE] # Left-to-right, TE mode.
		eigenvec2 = [ExTE, EyTE, EzTEm, HxTEm, HyTEm, HzTE] # Right-to-left, TE mode.
		eigenvec3 = [ExTM, EyTM, EzTMp, HxTMp, HyTMp, HzTM] # Left-to-right, TM mode.
		eigenvec4 = [ExTM, EyTM, EzTMm, HxTMm, HyTMm, HzTM] # Right-to-left, TM mode.

		eigvalues = [kzTEp, kzTEm, kzTMp, kzTMm]
		eigenvecs = hcat(eigenvec1, eigenvec2, eigenvec3, eigenvec4)

		return eigvalues, eigenvecs
	end

	"""
		_make_WhVh(ω, kx, ky, kz, eigvecs, z)

	Make Wh and Vh matrices.

    # Arguments
	- 'ω::Real': angular frequency of the input wave.
	- 'kx::Real': kx in free space.
	- 'ky::Real': ky in free space.
	- 'kz::Real': kz in free space.

    # Returns
	- 'Wh::AbstractMatrix': 2x2 identity matrix.
	- 'Vh::AbstractMatrix': 2x2 matrix.

    # Examples
	'''julia
	julia>
	'''

    # Notes
	"""
	function _make_WhVh(ω::Real, kx::Real, ky::Real, kz::Real)

		μ_0 = 4*π*10^-7

		Wh = la.I # The most Julianic way of expressing the identity matrix.
		Vh = zeros(ComplexF64, 2, 2)

		Vh[1,1] = kx * ky / kz / ω / μ_0
		Vh[1,2] = (kx^2 + kz^2) / kz / ω / μ_0
		Vh[2,1] =-(ky^2 + kz^2) / kz / ω / μ_0
		Vh[2,2] =-kx * ky / kz / ω / μ_0

		return Wh, Vh
	end

	"""
		_make_WpVp(eigvecs, z)

	Make some matrices and calculate the coupling coeffient matrix operators.

    # Arguments

    # Returns

    # Examples
	```julia
	julia>
	```

    # Notes
	"""
	function _make_WpVp(
		eigvalues::AbstractVector, 
		eigvecs::AbstractMatrix, 
		z::Real
	)

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

	function _make_WmVm(
		eigvalues::AbstractVector, 
		eigvecs::AbstractMatrix, 
		z::Real
	)

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

	function get_the_left_to_right_operators(
		ω::Number, 
		kx0::Number, 
		ky0::Number, 
		kz0::Number, 
		eigvalues::AbstractVector, 
		eigvectors::AbstractMatrix, 
		zm::Real, 
		zp::Real
		)
		
		Wh, Vh = _make_WhVh(ω, kx0, ky0, kz0)
		Wp0, Vp0 = _make_WpVp(eigvalues, eigvectors, 0.)
		Wpp, Vpp = _make_WpVp(eigvalues, eigvectors, zp-zm)
		Wm0, Vm0 = _make_WmVm(eigvalues, eigvectors, 0)
		Wmm, Vmm = _make_WmVm(eigvalues, eigvectors, zm-zp)

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

		R = inv(Wh) * (Wp0*cap + Wmm*cam - Wh*la.I)
		T = inv(Wh) * (Wpp*cap + Wm0*cam)

		return cap, cam, R, T
	end

	function get_the_right_to_left_operators(
		ω::Number, 
		kx0::Number, 
		ky0::Number, 
		kz0::Number, 
		eigvalues::AbstractVector, 
		eigvectors::AbstractMatrix, 
		zm::Real, 
		zp::Real
		)
		
		Wh, Vh = _make_WhVh(ω, kx0, ky0, kz0)
		Wp0, Vp0 = _make_WpVp(eigvalues, eigvectors, 0.)
		Wpp, Vpp = _make_WpVp(eigvalues, eigvectors, zp-zm)
		Wm0, Vm0 = _make_WmVm(eigvalues, eigvectors, 0)
		Wmm, Vmm = _make_WmVm(eigvalues, eigvectors, zm-zp)

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

		R = inv(Wh) * (Wp0*cbp + Wm0*cbm - Wh*la.I)
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

		kzns, eigvectors = TMM.get_eigenvectors(k0, kx0, ky0, mur, epsr)

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

		cap, cam, R1, T1 = TMM.get_the_left_to_right_operators(ω, kx0, ky0, kz0, kzns, eigvectors, zm, zp)

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

		kzns, eigvectors = TMM.get_eigenvectors(k0, kx0, ky0, mur, epsr)

		kzTEp = kzns[1]
		kzTEm = kzns[2]
		kzTMp = kzns[3]
		kzTMm = kzns[4]

		eigTEp = eigvectors[:,1]
		eigTEm = eigvectors[:,2]
		eigTMp = eigvectors[:,3]
		eigTMm = eigvectors[:,4]

		cbp, cbm, R, T = TMM.get_the_right_to_left_operators(ω, kx0, ky0, kz0, kzns, eigvectors, zm, zp)

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

	"""
		get_scc_of_a_block(ω, k, θ, ϕ, mur, epsr, zm, zp)

	Obtain the scattering matrix and the coupling coefficient matrix operators
	of a single block that is surrounded by the free space.
	
	# Arguments

	# Returns

	# Examples
	"""
	function get_scc_of_a_block(ω, k0, θ, ϕ, mur, epsr, zm, zp)

		# wave vector in a medium.
		kx0 = k0 * sin(θ) * cos(ϕ)
		ky0 = k0 * sin(θ) * sin(ϕ)
		kz0 = k0 * cos(θ)

		eigvals, eigvecs = TMM.get_eigenvectors(k0, kx0, ky0, mur, epsr)
		cap, cam, Ra, Ta = TMM.get_the_left_to_right_operators(ω, kx0, ky0, kz0, eigvals, eigvecs, zm, zp)
		cbp, cbm, Rb, Tb = TMM.get_the_right_to_left_operators(ω, kx0, ky0, kz0, eigvals, eigvecs, zm, zp)

		S = zeros(ComplexF64, 4, 4)
		S[1:2, 1:2] = Ta
		S[1:2, 3:4] = Ra
		S[3:4, 1:2] = Rb
		S[3:4, 3:4] = Tb

		Ca = zeros(ComplexF64, 4, 2)
		Cb = zeros(ComplexF64, 4, 2)

		Ca[1:2, :] = cap
		Ca[3:4, :] = cam
		Cb[1:2, :] = cbp
		Cb[3:4, :] = cbm

		return S, Ca, Cb

	end

	function get_scc_of_each_n_blocks(N, ω, k0, θ, ϕ, zs, murs, epsrs)

		@assert length(zs) == (length(murs)+1) == (length(epsrs)+1) "Please insert the "

		# For each n blocks, get the S matrix and
		# the Coupling coefficient matrices.
		Ss = Matrix{ComplexF64}[]
		Cas = Matrix{ComplexF64}[]
		Cbs = Matrix{ComplexF64}[]
		# Ras = Matrix{Float64}[]
		# Rbs = Matrix{Float64}[]
		# Tas = Matrix{Float64}[]
		# Tbs = Matrix{Float64}[]

		for i in 1:N

			S, Ca, Cb = get_scc_of_a_block(ω, k0, θ, ϕ, murs[i], epsrs[i], zs[i], zs[i+1])

			push!(Ss, S)
			push!(Cas, Ca)
			push!(Cbs, Cb)
			# push!(Tas, Ta)
			# push!(Ras, Ra)
			# push!(Rbs, Rb)
			# push!(Tbs, Tb)
		end

		return Ss, Cas, Cbs

	end

	function redheffer(S1, S2)

		T11_LtoR = S1[1:2, 1:2]
		R11_LtoR = S1[1:2, 3:4]
		R11_RtoL = S1[3:4, 1:2]
		T11_RtoL = S1[3:4, 3:4]
		
		T22_LtoR = S2[1:2, 1:2]
		R22_LtoR = S2[1:2, 3:4]
		R22_RtoL = S2[3:4, 1:2]
		T22_RtoL = S2[3:4, 3:4]

		R12_RtoL = R11_RtoL + T11_RtoL * inv(la.I - R22_RtoL*R11_LtoR) * R22_RtoL * T11_LtoR
		T12_LtoR = T22_LtoR * inv(la.I - R11_LtoR * R22_RtoL) * T11_LtoR
		R12_LtoR = R22_LtoR + T22_LtoR * inv(la.I - R11_LtoR*R22_RtoL) * R11_LtoR * T22_RtoL
		T12_RtoL = T11_RtoL * inv(la.I - R22_RtoL * R11_LtoR) * T22_RtoL

		S12 = zeros(ComplexF64, 4, 4)
		S12[1:2, 1:2] = T12_LtoR
		S12[1:2, 3:4] = R12_LtoR
		S12[3:4, 1:2] = R12_RtoL
		S12[3:4, 3:4] = T12_RtoL

		return S12

	end

	function redheffer_two_blocks(S1, S2, Ca1, Cb1, Ca2, Cb2)

		S12 = redheffer(S1, S2)

		T11_LtoR = S1[1:2, 1:2]
		R11_LtoR = S1[1:2, 3:4]
		R11_RtoL = S1[3:4, 1:2]
		T11_RtoL = S1[3:4, 3:4]
		
		T22_LtoR = S2[1:2, 1:2]
		R22_LtoR = S2[1:2, 3:4]
		R22_RtoL = S2[3:4, 1:2]
		T22_RtoL = S2[3:4, 3:4]

		cap111 = Ca1[1:2, :]
		cam111 = Ca1[3:4, :]
		cbp111 = Cb1[1:2, :]
		cbm111 = Cb1[3:4, :]

		cap222 = Ca2[1:2, :]
		cam222 = Ca2[3:4, :]
		cbp222 = Cb2[1:2, :]
		cbm222 = Cb2[3:4, :]

		cap121 = cap111 + cbp111 * inv(la.I-R22_RtoL*R11_LtoR) * R22_RtoL * T11_LtoR
		cam121 = cam111 + cbm111 * inv(la.I-R22_RtoL*R11_LtoR) * R22_RtoL * T11_LtoR

		cbp121 = cbp111 * inv(la.I - R22_RtoL * R11_LtoR) * T22_RtoL
		cbm121 = cbm111 * inv(la.I - R22_RtoL * R11_LtoR) * T22_RtoL

		cap122 = cap222 * inv(la.I - R11_LtoR*R22_RtoL) * T11_LtoR
		cam122 = cam222 * inv(la.I - R11_LtoR*R22_RtoL) * T11_LtoR

		cbp122 = cbp222 + cap222 * inv(la.I - R11_LtoR*R22_RtoL) * R11_LtoR * T22_RtoL
		cbm122 = cbm222 + cam222 * inv(la.I - R11_LtoR*R22_RtoL) * R11_LtoR * T22_RtoL

		ca121 = zeros(ComplexF64, 4, 2)
		ca122 = zeros(ComplexF64, 4, 2)
		cb121 = zeros(ComplexF64, 4, 2)
		cb122 = zeros(ComplexF64, 4, 2)

		ca121[1:2, :] = cap121
		ca121[3:4, :] = cam121
		
		cb121[1:2, :] = cbp121
		cb121[3:4, :] = cbm121
		
		ca122[1:2, :] = cap122
		ca122[3:4, :] = cam122
		
		cb122[1:2, :] = cbp122
		cb122[3:4, :] = cbm122

		Ca = [ca121, ca122] # of which the type is 'Vector{Matrix{ComplexF64}}'.
		Cb = [cb121, cb122]
		
		return S12, Ca, Cb
	end

	"""
		redheffer_left_blocks(m, l, Sl, Sr, Cals, Cbls, Cars, Cbrs)

	Obtain the scattering matrix and the coupling coefficient matrices
	of the new blocks using the previous partial blocks.
	The previous partial blocks on the left will be denoted as 'l' 
	in the name of the variable and the one the right will be 'r'.

	# Arguments
	- 'm':
	- 'l':
	- 'Sl':
	- 'Sr':
	- 'Cals':
	- 'Cbls':
	- 'Cars':
	- 'Cbrs':

	# Returns

	# Examples

	"""
	function redheffer_L_blocks(
		m::Int,
		Sl::Matrix{ComplexF64}, 
		Sr::Matrix{ComplexF64}, 
		Cal::Matrix{ComplexF64},
		Cbl::Matrix{ComplexF64}
	)

		Tl_LtoR = Sl[1:2, 1:2]
		Rl_LtoR = Sl[1:2, 3:4]
		Rl_RtoL = Sl[3:4, 1:2]
		Tl_RtoL = Sl[3:4, 3:4]

		Tr_LtoR = Sr[1:2, 1:2]
		Rr_LtoR = Sr[1:2, 3:4]
		Rr_RtoL = Sr[3:4, 1:2]
		Tr_RtoL = Sr[3:4, 3:4]

		# Obtain the reflection and transmission coefficient matrix operators.
		R_RtoL_nnml = Rl_RtoL + Tl_RtoL*inv(la.I - Rr_RtoL*Rl_LtoR)*Rr_RtoL*Tl_LtoR
		T_LtoR_nnml = Tr_LtoR * inv(la.I - Rl_LtoR*Rr_RtoL)*Tl_LtoR

		for m in 1:(m+1)
			# Obtain the coupling coefficient matrix operators of each block.
			ca_nnml_L = zeros(ComplexF64, 4, 2)
			cb_nnml_L = zeros(ComplexF64, 4, 2)

			calkp = Cal[1:2, :]
			calkm = Cal[3:4, :]
			cblkp = Cbl[1:2, :]
			cblkm = Cbl[3:4, :]

			cap = calkp + cblkp*inv(la.I - Rr_RtoL*Rl_LtoR)*Rr_RtoL*Tl_LtoR
			cam = calkm + cblkm*inv(la.I - Rr_RtoL*Rl_LtoR)*Rr_RtoL*Tl_LtoR
			cbp = cblkp * inv(la.I - Rr_RtoL*Rl_LtoR)*Tr_RtoL
			cbm = cblkm * inv(la.I - Rr_RtoL*Rl_LtoR)*Tr_RtoL

			ca_nnml_L[1:2, :] = cap
			ca_nnml_L[3:4, :] = cam
			cb_nnml_L[1:2, :] = cbp
			cb_nnml_L[3:4, :] = cbm
		end

		return R_RtoL_nnml, T_LtoR_nnml, ca_nnml_L, cb_nnml_L
	end

	function redheffer_R_blocks(
		l::Int,
		Sl::Matrix{ComplexF64}, 
		Sr::Matrix{ComplexF64}, 
		Car::Matrix{ComplexF64}, 
		Cbr::Matrix{ComplexF64}
	)

		Tl_LtoR = Sl[1:2, 1:2]
		Rl_LtoR = Sl[1:2, 3:4]
		Rl_RtoL = Sl[3:4, 1:2]
		Tl_RtoL = Sl[3:4, 3:4]

		Tr_LtoR = Sr[1:2, 1:2]
		Rr_LtoR = Sr[1:2, 3:4]
		Rr_RtoL = Sr[3:4, 1:2]
		Tr_RtoL = Sr[3:4, 3:4]

		# Obtain the reflection and transmission coefficient matrix operators.
		R_LtoR_nnml = Rr_LtoR + Tr_LtoR*inv(la.I - Rl_LtoR*Rr_RtoL)*Rl_LtoR*Tr_RtoL
		T_RtoL_nnml = Tl_RtoL * inv(la.I - Rr_RtoL*Rl_LtoR)*Tr_RtoL

		for l in 1:l
			# Obtain the coupling coefficient matrix operators of each block.
			ca_nnml_R = zeros(ComplexF64, 4, 2)
			cb_nnml_R = zeros(ComplexF64, 4, 2)

			carkp = Car[1:2, :]
			carkm = Car[3:4, :]
			cbrkp = Cbr[1:2, :]
			cbrkm = Cbr[3:4, :]

			cap = carkp * inv(la.I - Rl_LtoR*Rr_RtoL) * Tl_LtoR
			cam = carkm * inv(la.I - Rl_LtoR*Rr_RtoL) * Tl_LtoR
			cbp = cbrkp + carkp * inv(la.I-Rl_LtoR*Rr_RtoL)*Rl_LtoR*Tr_RtoL
			cbm = cbrkm + carkm * inv(la.I-Rl_LtoR*Rr_RtoL)*Rl_LtoR*Tr_RtoL

			ca_nnml_R[1:2, :] = cap
			ca_nnml_R[3:4, :] = cam
			cb_nnml_R[1:2, :] = cbp
			cb_nnml_R[3:4, :] = cbm

		end

		return R_LtoR_nnml, T_RtoL_nnml, ca_nnml_R, cb_nnml_R
	end

	"""
		redheffer_n_blocks(Ss, Cas, Cbs)

	Obtain the scattering matrix and the coupling coefficient matrix
	operators of the combined multiblock M^(n, n+m+l).

	# Arguments
	# Returns
	# Examples
	"""
	function redheffer_n_blocks(
		ω::Real, 
		θ::Real, 
		ϕ::Real,
		zs::Vector,
		murs::Vector, 
		epsrs::Vector
	)

		c = 299792458 # m/s.
		freq = ω / 2 / π
		λ0 = c / freq # wavelength in free space.
		k0 = 2 * π / λ0

		# The total number of the blocks.
		N = length(murs) # It should be noted that N = m + l + 1
		@assert N == length(epsrs) "The length of mur and epsr are not the same."

		Ss, Cas, Cbs = get_scc_of_each_n_blocks(N, ω, k0, θ, ϕ, zs, murs, epsrs)

		if N == 1
			return Ss[1], Cas[1], Cbs[1]
		elseif N == 2
			S12, Ca12, Cb12 = redheffer_two_blocks(Ss[1], Ss[2], Cas[1], Cbs[1], Cas[2], Cbs[2])
			return S12, Ca12, Cb12
		elseif N > 2
			# For more than one block, let us consider the left and right portion of the blocks.
			n = 1 # The index of the first slab on the left.
			m = round(Int, N/2-1) # the number of the blocks on the left - 1.
			l = N - m - 1 # the number of the blocks on the right.

			combinedSs = Vector{Matrix{ComplexF64}}()
			combinedCas = Vector{Vector{Matrix{ComplexF64}}}() 
			combinedCbs = Vector{Vector{Matrix{ComplexF64}}}() 

			S12, Ca12, Cb12 = redheffer_two_blocks(Ss[1], Ss[2], Cas[1], Cbs[1], Cas[2], Cbs[2])

			push!(combinedSs, S12)   # Ss = [S^(1,1), S^(1,2), S^(1,3), ...]
			push!(combinedCas, Ca12) # Cas = [Ca^(n,n+m+l)_(n,n), Ca^(n,n+m+l)_(1,2), ...]
			push!(combinedCbs, Cb12) # Cbs = [Cb^(n,n+m+l)_(n,n), Cb^(n,n+m+l)_(1,2), ...]

			for i in 3:N

				newS = zeros(ComplexF64, 4, 4)

				Tl_LtoR = combinedSs[end][1:2, 1:2]
				Rl_LtoR = combinedSs[end][1:2, 3:4]
				Rl_RtoL = combinedSs[end][3:4, 1:2]
				Tl_RtoL = combinedSs[end][3:4, 3:4]

				Tr_LtoR = Ss[i][1:2, 1:2]
				Rr_LtoR = Ss[i][1:2, 3:4]
				Rr_RtoL = Ss[i][3:4, 1:2]
				Tr_RtoL = Ss[i][3:4, 3:4]

				# Obtain the reflection and transmission coefficient matrix operators.
				T_LtoR_nnml = Tr_LtoR * inv(la.I - Rl_LtoR*Rr_RtoL)*Tl_LtoR
				T_RtoL_nnml = Tl_RtoL * inv(la.I - Rr_RtoL*Rl_LtoR)*Tr_RtoL
				R_LtoR_nnml = Rr_LtoR + Tr_LtoR*inv(la.I - Rl_LtoR*Rr_RtoL)*Rl_LtoR*Tr_RtoL
				R_RtoL_nnml = Rl_RtoL + Tl_RtoL*inv(la.I - Rr_RtoL*Rl_LtoR)*Rr_RtoL*Tl_LtoR

				# Update the new scattering matrix.
				newS[1:2, 1:2] = T_LtoR_nnml
				newS[1:2, 3:4] = R_LtoR_nnml
				newS[3:4, 1:2] = R_RtoL_nnml
				newS[3:4, 3:4] = T_RtoL_nnml

				push!(combinedSs, newS)

				# R_RtoL_nnml, T_LtoR_nnml, ca_nnmls_L, cb_nnmls_L = redheffer_L_blocks(m, combinedSs[end], Ss[i+1], combinedCas[end], Cas[i+1])
				# R_LtoR_nnml, T_RtoL_nnml, ca_nnmls_R, cb_nnmls_R = redheffer_R_blocks(l, combinedSs[end], Ss[i+1], combinedCbs[end], Cbs[i+1])

				ca_nnml_LRs = Vector{Matrix{ComplexF64}}(undef, i)
				cb_nnml_LRs = Vector{Matrix{ComplexF64}}(undef, i)

				# Obtain the coupling coefficient matrices for the left blocks.
				for k in 1:(i-1)

					last_combinedCa = combinedCas[end]
					last_combinedCb = combinedCbs[end]

					calp = last_combinedCa[k][1:2, :]
					calm = last_combinedCa[k][3:4, :]
					cblp = last_combinedCb[k][1:2, :]
					cblm = last_combinedCb[k][3:4, :]

					cap = calp + cblp*inv(la.I - Rr_RtoL*Rl_LtoR)*Rr_RtoL*Tl_LtoR
					cam = calm + cblm*inv(la.I - Rr_RtoL*Rl_LtoR)*Rr_RtoL*Tl_LtoR
					cbp = cblp * inv(la.I - Rr_RtoL*Rl_LtoR)*Tr_RtoL
					cbm = cblm * inv(la.I - Rr_RtoL*Rl_LtoR)*Tr_RtoL

					newly_combinedCa_Ls = zeros(ComplexF64, 4, 2)
					newly_combinedCb_Ls = zeros(ComplexF64, 4, 2)

					newly_combinedCa_Ls[1:2, :] = cap
					newly_combinedCa_Ls[3:4, :] = cam
					newly_combinedCb_Ls[1:2, :] = cbp
					newly_combinedCb_Ls[3:4, :] = cbm

					ca_nnml_LRs[k] = newly_combinedCa_Ls
					cb_nnml_LRs[k] = newly_combinedCb_Ls
				end

				# Obtain the coupling coefficient matrix for a right block.
				newly_combinedCa_R = zeros(ComplexF64, 4, 2)
				newly_combinedCb_R = zeros(ComplexF64, 4, 2)

				ca_to_combine = Cas[i]
				cb_to_combine = Cbs[i]

				carp = ca_to_combine[1:2, :]
				carm = ca_to_combine[3:4, :]
				cbrp = cb_to_combine[1:2, :]
				cbrm = cb_to_combine[3:4, :]

				cap = carp * inv(la.I - Rl_LtoR*Rr_RtoL) * Tl_LtoR
				cam = carm * inv(la.I - Rl_LtoR*Rr_RtoL) * Tl_LtoR
				cbp = cbrp + carp * inv(la.I-Rl_LtoR*Rr_RtoL)*Rl_LtoR*Tr_RtoL
				cbm = cbrm + carm * inv(la.I-Rl_LtoR*Rr_RtoL)*Rl_LtoR*Tr_RtoL

				newly_combinedCa_R[1:2, :] = cap
				newly_combinedCa_R[3:4, :] = cam
				newly_combinedCb_R[1:2, :] = cbp
				newly_combinedCb_R[3:4, :] = cbm

				ca_nnml_LRs[i] = newly_combinedCa_R
				cb_nnml_LRs[i] = newly_combinedCb_R

				push!(combinedCas, ca_nnml_LRs)
				push!(combinedCbs, cb_nnml_LRs)
			end

			# for i in (m+1):(n+m+l)
			# end

			return combinedSs[end], combinedCas[end], combinedCbs[end]
		else
			throw(DomainError(N, "negative N not allowed."))
		end

	end

	"""
		redheffer_left_half_infinite()
	
	Perform redheffer star product to obtain the scattering matrix and
	the coupling coefficient matrix operators on the left half-infinite
	block.

	# Arguments
	# Returns
	# Examples
	"""
	function redheffer_left_half_infinite(ω, θ, ϕ, zm)

		mur_lhi = 1
		epsr_lhi = 1
		n_lhi = sqrt(mur_lhi*epsr_lhi) # Note that we assumed the surrounding medium is vaccuum of air.
		c = 299792458 # m/s.
		freq = ω / 2 / π
		λn = c / n_lhi / freq
		k0 = 2 * π / λn

		# wave vector in free space.
		kx0 = k0 * sin(θ) * cos(ϕ)
		ky0 = k0 * sin(θ) * sin(ϕ)
		kz0 = k0 * cos(θ)

		eigvals0, eigvecs0 = TMM.get_eigenvectors(k0, kx0, ky0, mur_lhi, epsr_lhi)
		Wh, Vh = TMM._make_WhVh(ω, kx0, ky0, kz0)
		Wp, Vp = TMM._make_WpVp(eigvals0, eigvecs0, zm)
		Wm, Vm = TMM._make_WmVm(eigvals0, eigvecs0, zm)

		A = inv(Wh)*Wm - inv(Vh)*Vm
		B = inv(Wh)*Wp - inv(Vh)*Vp
		C = inv(Wm)*Wh - inv(Vm)*Vh
		D = inv(Wm)*Wp - inv(Vm)*Vp
		E = inv(Wm)*Wh + inv(Vm)*Vh

		R_RtoL =  -inv(A) * B 
		T_LtoR =   inv(C) * D
		T_RtoL = 2*inv(A)
		R_LtoR =  -inv(C) * E

		S = zeros(ComplexF64, 4, 4)
		S[1:2, 1:2] = T_LtoR
		S[1:2, 3:4] = R_LtoR
		S[3:4, 1:2] = R_RtoL
		S[3:4, 3:4] = T_RtoL

		return S

	end

	function redheffer_right_half_infinite(ω, θ, ϕ, zp)
		
		n_rhi = 1 # Note that we assumed the surrounding medium is vaccuum of air.
		c = 299792458 # m/s.
		freq = ω / 2 / π
		λn = c / n_rhi / freq
		k0 = 2 * π / λn

		# wave vector in free space.
		kx0 = k0 * sin(θ) * cos(ϕ)
		ky0 = k0 * sin(θ) * sin(ϕ)
		kz0 = k0 * cos(θ)

		eigvals0, eigvecs0 = TMM.get_eigenvectors(k0, kx0, ky0, 1, 1)
		Wh, Vh = TMM._make_WhVh(ω, kx0, ky0, kz0)
		Wp, Vp = TMM._make_WpVp(eigvals0, eigvecs0, zp)
		Wm, Vm = TMM._make_WmVm(eigvals0, eigvecs0, zp)

		A = inv(Wh)*Wp + inv(Vh)*Vp
		B = inv(Wh)*Wh + inv(Vp)*Vh
		C = inv(Wp)*Wh - inv(Vp)*Vh
		D = inv(Wp)*Wh + inv(Vp)*Vh
		E = inv(Wp)*Wm - inv(Vp)*Vm
		F = inv(Wh)*Wm + inv(Vh)*Vm

		T_LtoR = 2*inv(A)
		R_RtoL =  -inv(B) * C 
		T_RtoL =   inv(D) * E
		R_LtoR =  -inv(A) * F

		S = zeros(ComplexF64, 4, 4)
		S[1:2, 1:2] = T_LtoR
		S[1:2, 3:4] = R_LtoR
		S[3:4, 1:2] = R_RtoL
		S[3:4, 3:4] = T_RtoL

		return S

	end

	function redheffer_combine_all_blocks(ω, θ, ϕ, zs, murs, epsrs)

		zm = zs[1]
		zp = zs[end]

		Slhi = redheffer_left_half_infinite(ω, θ, ϕ, zm)
		Srhi = redheffer_right_half_infinite(ω, θ, ϕ, zp)

		T00_LtoR = Slhi[1:2, 1:2]
		R00_LtoR = Slhi[1:2, 3:4]
		R00_RtoL = Slhi[3:4, 1:2]
		T00_RtoL = Slhi[3:4, 3:4]

		Tn1n1_LtoR = Srhi[1:2, 1:2]
		Rn1n1_LtoR = Srhi[1:2, 3:4]
		Rn1n1_RtoL = Srhi[3:4, 1:2]
		Tn1n1_RtoL = Srhi[3:4, 3:4]

		S1n, Ca1n, Cb1n = redheffer_n_blocks(ω, θ, ϕ, zs, murs, epsrs)

		# println(size(S1n))
		# println(size(Ca1n))
		# println(size(Cb1n))
		# println(size(Ca1n[1]))
		# println(size(Cb1n[2]))

		T1n_LtoR = S1n[1:2, 1:2]
		R1n_LtoR = S1n[1:2, 3:4]
		R1n_RtoL = S1n[3:4, 1:2]
		T1n_RtoL = S1n[3:4, 3:4]

		S0n = redheffer(Slhi, S1n)
		S0n1 = redheffer(S0n, Srhi)

		# println(size(S0n1))

		T0n_LtoR = S0n[1:2, 1:2]
		R0n_LtoR = S0n[1:2, 3:4]
		R0n_RtoL = S0n[3:4, 1:2]
		T0n_RtoL = S0n[3:4, 3:4]

		T0n1_LtoR = S0n1[1:2, 1:2]
		R0n1_LtoR = S0n1[1:2, 3:4]
		R0n1_RtoL = S0n1[3:4, 1:2]
		T0n1_RtoL = S0n1[3:4, 3:4]

		# Combine the left half infinite block.
		N = length(Ca1n)
		Ca0n = Vector{Matrix{ComplexF64}}(undef, N)
		Cb0n = Vector{Matrix{ComplexF64}}(undef, N)
		Ca0n1 = Vector{Matrix{ComplexF64}}(undef, N)
		Cb0n1 = Vector{Matrix{ComplexF64}}(undef, N)

		println(size(Ca0n))
		println(typeof(Ca0n))

		for k in 1:N
			# combine the left half infinite block.
			Ca0n[k] = Ca1n[k] * (inv(la.I - R00_LtoR*R1n_RtoL) * T00_LtoR)
			Cb0n[k] = Cb1n[k] + Ca1n[k] * inv(la.I - R00_LtoR*R1n_RtoL) * R00_LtoR * T1n_RtoL

			# combine the right half infinite block.
			Ca0n1[k] = Ca0n[k] + (Cb0n[k] * inv(la.I - Rn1n1_RtoL*R0n_LtoR) * Rn1n1_RtoL * T0n_LtoR)
			Cb0n1[k] = Cb0n[k] * inv(la.I - Rn1n1_RtoL*R0n_LtoR) * Tn1n1_RtoL

		end

		return S0n1, Ca0n1, Cb0n1

	end

end # end of a module.