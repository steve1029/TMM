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
	export get_scc_of_a_block
	export get_scc_of_each_n_blocks
	export redheffer_two_blocks
	export redheffer_L_blocks
	export redheffer_R_blocks
	export redheffer_n_blocks

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
	function get_eigenvectors(k0, kx0, ky0, mur, epsr, impedance)

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
	function _make_WhVh(μ_0, ω::Real, kx::Real, ky::Real, kz::Real)

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
		μ_0::Real,
		ω::Number, 
		kx0::Number, 
		ky0::Number, 
		kz0::Number, 
		eigvalues::AbstractVector, 
		eigvectors::AbstractMatrix, 
		zm::Real, 
		zp::Real
		)
		
		Wh, Vh = _make_WhVh(μ_0, ω, kx0, ky0, kz0)
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
		μ_0::Real,
		ω::Number, 
		kx0::Number, 
		ky0::Number, 
		kz0::Number, 
		eigvalues::AbstractVector, 
		eigvectors::AbstractMatrix, 
		zm::Real, 
		zp::Real
		)
		
		Wh, Vh = _make_WhVh(μ_0, ω, kx0, ky0, kz0)
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

	function get_scc_of_a_block(μ_0, ω, k0, θ, ϕ, mur, epsr, impedance, zm, zp)

		# wave vector in free space.
		kx0 = k0 * sin(θ) * cos(ϕ)
		ky0 = k0 * sin(θ) * sin(ϕ)
		kz0 = k0 * cos(θ)

		eigvals, eigvecs = TMM.get_eigenvectors(k0, kx0, ky0, mur, epsr, impedance)
		cap, cam, Ra, Ta = TMM.get_the_left_to_right_operators(μ_0, ω, kx0, ky0, kz0, eigvals, eigvecs, zm, zp)
		cbp, cbm, Rb, Tb = TMM.get_the_right_to_left_operators(μ_0, ω, kx0, ky0, kz0, eigvals, eigvecs, zm, zp)

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

	function get_scc_of_each_n_blocks(N, μ_0, ω, impedance, k0, θ, ϕ, zs, murs, epsrs)

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

			S, Ca, Cb = get_scc_of_a_block(μ_0, ω, k0, θ, ϕ, murs[i], epsrs[i], impedance, zs[i], zs[i+1])

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

	function redheffer_two_blocks(S1, S2, Ca1, Cb1, Ca2, Cb2)

		T11a = S1[1:2, 1:2]
		R11a = S1[1:2, 3:4]
		T11b = S1[1:2, 1:2]
		R11b = S1[1:2, 3:4]
		
		T22a = S2[1:2, 1:2]
		R22a = S2[1:2, 3:4]
		T22b = S2[1:2, 1:2]
		R22b = S2[1:2, 3:4]

		R12a = R22a + T22a * inv(la.I - R11a*R22b) * R11a * T22b
		R12b = R11b + T11b * inv(la.I - R22b*R11a) * R22b * T11a
		T12a = T22a * inv(la.I - R11a * R22b) * T11a
		T12b = T11b * inv(la.I - R22b * R11a) * T22b

		S12 = zeros(ComplexF64, 4, 4)
		S12[1:2, 1:2] = T12a
		S12[1:2, 3:4] = R12a
		S12[3:4, 1:2] = T12b
		S12[3:4, 3:4] = R12b

		cap111 = Ca1[1:2, :]
		cam111 = Ca1[3:4, :]
		cbp111 = Cb1[1:2, :]
		cbm111 = Cb1[3:4, :]

		cap222 = Ca2[1:2, :]
		cam222 = Ca2[3:4, :]
		cbp222 = Cb2[1:2, :]
		cbm222 = Cb2[3:4, :]

		cap121 = cap111 + cbp111 * inv(la.I-R22b*R11a) * R22b * T11a
		cam121 = cam111 + cbm111 * inv(la.I-R22b*R11a) * R22b * T11a

		cbp121 = cbp111 * inv(la.I - R22b * R11a) * T22b
		cbm121 = cbm111 * inv(la.I - R22b * R11a) * T22b

		cap122 = cap222 * inv(la.I - R11a*R22b) * T11a
		cam122 = cam222 * inv(la.I - R11a*R22b) * T11a

		cbp122 = cbp222 + cap222 * inv(la.I - R11a*R22b) * R11a * T22b
		cbm122 = cbm222 + cam222 * inv(la.I - R11a*R22b) * R11a * T22b

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
		N::Int,
		μ_0::Real, 
		ω::Real, 
		impedance::Real,
		k0::Real,
		θ::Real, 
		ϕ::Real,
		zs::Vector,
		murs::Vector, 
		epsrs::Vector
	)

		# wave vector in free space.
		kx0 = k0 * sin(θ) * cos(ϕ)
		ky0 = k0 * sin(θ) * sin(ϕ)
		kz0 = k0 * cos(θ)

		Ss, Cas, Cbs = get_scc_of_each_n_blocks(N, μ_0, ω, impedance, k0, θ, ϕ, zs, murs, epsrs)

		# The total number of the blocks.
		N = length(Ss) # It should be noted that N = m + l + 1

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

end # end of a module.