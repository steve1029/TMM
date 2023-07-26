module single_block

	# The following packages are required in this module.
	using LinearAlgebra

	export get_Amatrix, get_Bmatrix, 
	get_HxHy, get_Ez, get_Hz, get_eigenvectors, 
	make_WhVh, make_WpVp, make_WmVm,
	get_the_left_to_right_operators,
	get_the_right_to_left_operators

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

	function get_HxHy(k0::Number, kzn::Number, A::AbstractMatrix, input::AbstractVector, impedance)

		H11 = (1im .* kzn./k0) .* (inv(A) * input)

		Hy = H11[1]
		Hx = H11[2]

		Hx *= (1im / impedance)
		Hy *= (1im / impedance)

		return Hx, Hy
	end

	function get_Ez(ky_bar, kx_bar, Hy, Hx, epsr, impedence)
		Hy *= (impedence / 1im)
		Hx *= (impedence / 1im)
		Ez = 1im.*(ky_bar.*Hx .- kx_bar.*Hy) ./ epsr
		return Ez
	end

	function get_Hz(ky_bar, kx_bar, Ey, Ex, mur, impedance)
		Hz = 1im.*(ky_bar.*Ex .- kx_bar.*Ey) ./ mur
		Hz*= (1im / impedance)
		return Hz
	end

    # Eigen vector and eigen value calculation in a single block.
	function get_eigenvectors(k0, kxn, kyn, mur, epsr, impedance)

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

		eigvals, EyEx_eigvecs = LinearAlgebra.eigen(C)

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
		HxTEp, HyTEp = get_HxHy(k0, kzTEp, A, [EyTE, ExTE], impedance)
		HxTEm, HyTEm = get_HxHy(k0, kzTEm, A, [EyTE, ExTE], impedance)
		HxTMp, HyTMp = get_HxHy(k0, kzTMp, A, [EyTM, ExTM], impedance)
		HxTMm, HyTMm = get_HxHy(k0, kzTMm, A, [EyTM, ExTM], impedance)

		EzTEp = get_Ez(ky_bar, kx_bar, HyTEp, HxTEp, epsrz, impedance)
		EzTEm = get_Ez(ky_bar, kx_bar, HyTEm, HxTEm, epsrz, impedance)
		EzTMp = get_Ez(ky_bar, kx_bar, HyTMp, HxTMp, epsrz, impedance)
		EzTMm = get_Ez(ky_bar, kx_bar, HyTMm, HxTMm, epsrz, impedance)

		HzTE = get_Hz(ky_bar, kx_bar, EyTE, ExTE, murz, impedance)
		HzTM = get_Hz(ky_bar, kx_bar, EyTM, ExTM, murz, impedance)

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

		Wh = LinearAlgebra.I # The most Julianic way of expressing the identity matrix.
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

end # end of a module.
