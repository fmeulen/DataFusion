# Inference algorithm, define a simple transition kernel
function customkernel(θ, s::Symbol, scale)
	θ° = deepcopy(θ)
	θ°[s] += 2.0*scale*(rand()-0.5)
	θ°
end

"""
	define a multicariate transition kernel
	X is a Rnaodm Variable defined in Distributions
	keys are the keys of the dictionaries
	WARNING: you need to be careful with the order of the dicitionary (usually
	 ordered with alphabetic order but not when there are greek letters)
"""
function joint_update(θ, X, keys)
	pos = true
	θ° = deepcopy(θ)
	new_values =  rand(X)
	i = 1
	for T in keys
		θ°[T] = θ°[T] + new_values[i]
		i += 1
	end
	θ°
end

function adaptive_step(count_adap, ϵ, dt, lb, ub, target, name_variable)
	if count_adap <= target
		arrow = "↓"
		ϵ = max(ϵ - dt, lb)
	else
		ϵ = min(ϵ + dt, ub)
		arrow = "↑"
	end
	println("adaptation $name_variable $arrow : new value : $ϵ, acc_rate : $(count_adap)")
	ϵ, 0
end


"""
 update covariance matrix and mean Q, μ, given the new
 set of parameters θ at iteration N
"""
function update_cov(Q, μ, θ, N)
	old_sum_sq = (N-1)/N*Q + μ*μ'
    μ = N/(N+1)*μ + θ/(N+1)
    new_sum_sq = old_sum_sq + (θ*θ')/N
    Q = new_sum_sq - (N+1)/N*(μ*μ')
	return Q, μ
end

#standard_guid_prop_time_transf(tt)

function simple_inference(AuxLaw, recording, impgrid, θinit;
				sc =0.1, ρ=0.5, num_steps=10^4,save_skip_path = 1,
				save_skip_para = 1, λ_cov = 0.05, λ_Id = 0.03)

	# initializations
	tts = impgrid#OBS.setup_time_grids(recording, dt)
	ρρ = [ρ for _ in tts]
	PP = build_guid_prop(AuxLaw, recording, tts)
	PP° = deepcopy(PP)

	y1 = rand(recording.x0_prior) # just returns the starting point
	XX, WW, Wnr = rand(PP, y1)
	XX°, WW° = trajectory(PP)

	ll = loglikhd(PP, XX)
	paths = []
	θ = θinit
	θθ = [Float64[θ[:ᾱ]], Float64[θ[:ω̄]], Float64[θ[:ϕ̄]],  Float64[θ[:δ]], Float64[θ[:σ̄]] ] # θ-vals to be saved
 	# counters
	imp_a_r = 0
	#param_a_r = 0
	param_joint = 0
	diag_joint = 0

	# adaptation
	μ = collect(values(θ))
	Q = Matrix{Float64}(I,length(θinit),length(θinit))
	β = 0.5

	# mcmc
	burn_in_adapt = 0


	Id = Diagonal(Q)
	X = MvNormal(λ_cov^2*Q)
	Y = MvNormal(λ_Id^2*Id)

	# MCMC
	for i in 1:num_steps
		# impute a path
		_, ll° = rand!(PP, XX°, WW°, WW, ρρ, Val(:ll), y1; Wnr=Wnr)
		if rand() < exp(ll°-ll)
			XX, WW, XX°, WW° = XX°, WW°, XX, WW
			ll = ll°
			imp_a_r += 1
		end

		# update parameters ᾱ, ω1, σ̄, δ1
		if rand() > β
			θ° = joint_update(θ, X, collect(keys(θ)))
			diag = false
		else
			diag = true
			θ° = joint_update(θ, Y, collect(keys(θ)))
		end
		GP.set_parameters!(PP°, θ°)
		recompute_guiding_term!(PP°)
		_, ll° = GP.solve_and_ll!(XX°, WW, PP°, y1)
		if randexp() > (ll - ll°) # prior ratio- diff𝓝(C_μ, C_σ, θ°[:C], θ[:C]) - diff𝓝(σy_μ, σy_σ, θ°[:σy], θ[:σy])) && all(collect(values(θ°)).>0.0) # uniform updates have no contribution to ll
			XX, PP, θ, XX°, PP°, θ° = XX°, PP°, θ°, XX, PP, θ
			ll = ll°
			diag == false ? param_joint += 1 : diag_joint += 1
		end
		if i % save_skip_para == 0  append!(θθ[1], [θ[:ᾱ]]), append!(θθ[2], [θ[:ω̄]]),
									append!(θθ[3], [θ[:ϕ̄]]), append!(θθ[4], [θ[:δ]]), append!(θθ[5], [θ[:σ̄]]) end
									# ordering of parameters is (ᾱ, ω̄, ϕ̄, δ, σ̄)


		# train covariance matrix
		Q, μ = update_cov(Q, μ, collect(values(θ)), i)
		# println(Q)
		# println(isposdef(Q))
		#display info
		if i % 500 == 0
			println("$i. ll=$ll,ᾱ=$(θ[:ᾱ]), ω̄=$(θ[:ω̄]), ϕ̄=$(θ[:ϕ̄]), δ=$(θ[:δ]), σ̄=$(θ[:σ̄])")
			λ_cov, param_joint  = adaptive_step(param_joint/500/(1-β), λ_cov, 0.01, 0.00005, 2.0, 0.4, :λ_cov)
			λ_Id, diag_joint  = adaptive_step(diag_joint/500/(β), λ_Id, 0.01, 0.0001, 2.0, 0.4, :λ_diag)
			Y = MvNormal(λ_Id^2*Diagonal(Q))
			X = MvNormal(λ_cov^2*Q)
		end
		if i % 1000 == 0
			β = max(0.1, β - 0.01)
			println("new beta $β")
		end
		i % save_skip_path == 0 && append!(paths, [deepcopy(XX)])
	end
	paths, θθ, Q
end
