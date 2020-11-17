# Model: d X_t = - α (X_t - μ(t)) d t + σ d W_t
# μ(t) = ∑ᵢ ξᵢ * ϕᵢ(t)
# extrinsic noise with variances governed by ψ

struct DF
    α::Float64  # mean reversion par
    ξ::Vector{Float64}  #  pars in periodic drift function
    σ2::Float64    # squared diffusivity
    ψ::Vector{Float64}  # vars on observation equation
    t::Vector{Float64}  # t_i i ∈ 1...n
    Δ::Vector{Float64}  # Δ_i = t_i-t_{i-1}
    typeobs::Vector
	J::Int64			# truncation level Fourier series (so 2J+1 basis functions)
	η::Float64
end

parameters(𝒫::DF) = (𝒫.α, 𝒫.ξ, 𝒫.σ2, 𝒫.ψ)

struct ObsGroup
	ind1::Vector{Int64} # indices in y where measurement device 1 is used
	ind2::Vector{Int64} # indices in y where measurement device 2 is used
	y1::Vector 			# y where measurement device 1 is used
	y2::Vector			# y where measurement device 2 is used

	function ObsGroup(𝒫, y)
		y1 = SArray{Tuple{1},Float64,1,1}[]
		y2 = SArray{Tuple{1},Float64,1,1}[]
		ind1 = Int64[]
		ind2 = Int64[]
		for i ∈ eachindex(y)
			if 𝒫.typeobs[i]=="obs1"
				push!(y1,y[i])
				push!(ind1,i)
			elseif 𝒫.typeobs[i]=="obs2"
				push!(y2,y[i])
				push!(ind2,i)
			else
				push!(y1,SMatrix{1,1}(y[i][1]))
				push!(y2,SMatrix{1,1}(y[i][2]))
				push!(ind1,i)
				push!(ind2,i)
			end
		end
		new(ind1,ind2,y1,y2)
	end


end


# Model equations
# x[k] = A(k-1) * x[k-1] + a[k-1] +  N(0,Q[k-1])
# y[k] = H(k) * x[k] + N(0,R(k))

# We assume typeobs="obs1" (only meas device 1 used), typeobs="obs2" (only meas device 2 used),
# typeobs="obs3" (both meas devices 1 and 2 used)


#μ(t,ξ) =          pdf(Beta(ξ[1],ξ[2]),mod(t,1.0))#
# function ϕ(x, J)
# 	out = [1.0]
# 	for j ∈ 1:J
# 		push!(out, j^(-2)*cos(2π*j*x))
# 		push!(out, j^(-2)*sin(2π*j*x))
# 	end
# 	out
# end
#ϕ0(x) = 0<x<1 ? 2x*(x<0.5) + 2(1-x)*(x>=0.5) : 0.0

"""
	ϕ0(x)
Returns quadratic B-spline function with support on [0,1], evaluated at x
"""
function ϕ0(x)
	x = 3.0x
	if 0.0<= x<= 1.0
		return 0.5x^2
	elseif 1.0 <= x <= 2.0
		return 0.75-(x-1.5)^2
	elseif 2.0 <= x <= 3.0
		return 0.5(3.0-x)^2
	else
		return 0.0
	end
end



"""
	ϕ(x, J)
Returns vector of basis-functions, each evaluated at x
"""
function ϕ(x, J)
	x = mod(x, 1.0)
	out = [1.0]
	for j ∈ 1:J
		for k ∈ 0:2^(j-1)-1
			push!(out, ϕ0(2^(j-1)*x-k)/j)
		end
	end
	out
end
nbasis(J) = 2^(J) # 2*J+1



μ(t,ξ,J) = dot(ξ, ϕ(t, J))

implicit = true
if implicit
    A(k,𝒫) = k==0 ?  SMatrix{1,1}([1.0]) : SMatrix{1,1}( [(1.0 + 𝒫.α * 𝒫.Δ[k])^(-1)] )
    a(k,𝒫) = k==0 ?  (@SVector [0.0]) :   (@SVector [A(k,𝒫)[1,1] * 𝒫.α * μ(𝒫.t[k+1],𝒫.ξ, 𝒫.J) * 𝒫.Δ[k]   ])
else
    A(k,𝒫) = k==0 ?  SMatrix{1,1}([1.0]) : SMatrix{1,1}( [1.0 - 𝒫.α * 𝒫.Δ[k]] )
    a(k,𝒫) = k==0 ?  (@SVector [0.0]) :   (@SVector [𝒫.α *  μ(𝒫.t[k+1],𝒫.ξ) * 𝒫.Δ[k]   ])
end

Q(k,𝒫) = k==0 ?  SMatrix{1,1}([0.0]) :   SMatrix{1,1}( [𝒫.σ2 * 𝒫.Δ[k]] )

H(k,𝒫) =  𝒫.typeobs[k]=="obs3" ? SMatrix{2,1}([1.0 1.0]) : SMatrix{1,1}([1.0])

function R(k,𝒫)
	if   𝒫.typeobs[k]=="obs1"
    	return  SMatrix{1,1}([𝒫.η * 𝒫.ψ[1] * 𝒫.ψ[2]])
    elseif    𝒫.typeobs[k]=="obs2"
        return  SMatrix{1,1}([𝒫.ψ[2]])
    else
        return  SMatrix{2,2}([𝒫.η * 𝒫.ψ[1] * 𝒫.ψ[2] 0.0; 0.0 𝒫.ψ[2]])
    end
end


## mcmc updates for pars

function SS(𝒫, x)
	S = 0.0
	for k ∈ 2:length(x)
		S += (x[k] - A(k-1,𝒫) * x[k-1] - a(k-1,𝒫))[1]^2/𝒫.Δ[k-1]
	end
	S
end

function update_σ2(𝒫, x, Aσ, Bσ)
	m = length(x)-1
	S = SS(𝒫, x)
	rand(InverseGamma(0.5m+ Aσ, 0.5S + Bσ))
end

logtargetψ1(ψ1, ψ2 ,n ,S, η) = -0.5n * log(ψ1) - 0.5S/(η * ψ1 * ψ2)  # needed for MH-update first element of ψ

function update_ψ(𝒫::DF, 𝒢::ObsGroup, x, acc, prior_ψ1, Aψ, Bψ, propσ)
	S1 = norm(ec1(𝒢.y1-x[𝒢.ind1]))^2
	S2 = norm(ec1(𝒢.y2-x[𝒢.ind2]))^2
	n1 = length(𝒢.y1)
	n2 = length(𝒢.y2)
	ψ2 = rand(InverseGamma(Aψ + 0.5*(n1 + n2), Bψ + 0.5*(S1/(𝒫.η * 𝒫.ψ[1]) + S2) ))
	ψ1 = 𝒫.ψ[1]
	ψ1ᵒ = ψ1/(ψ1 + (1-ψ1)* exp(propσ*randn()))
	Δll = logtargetψ1(ψ1ᵒ,ψ2,n1,S1,𝒫.η) - logtargetψ1(ψ1,ψ2,n1,S1,𝒫.η)
	log_jacob_term = log(ψ1ᵒ*(1.0-ψ1ᵒ)) - log(ψ1*(1.0-ψ1))
	Υ = Δll + log_jacob_term  + logpdf(prior_ψ1, ψ1ᵒ) - logpdf(prior_ψ1, ψ1)
	if log(rand()) < Υ
		ψ1 = ψ1ᵒ
		acc[2] += 1
	end
	return [ψ1, ψ2], acc
end


function update_α(𝒫, x, acc, prior_α, propσ)
	α = 𝒫.α
	αᵒ = α * exp(propσ*randn())
	𝒫ᵒ = DF(αᵒ, 𝒫.ξ, 𝒫.σ2, 𝒫.ψ, 𝒫.t, 𝒫.Δ, 𝒫.typeobs, 𝒫.J, 𝒫.η)
	Δll = -0.5*(SS(𝒫ᵒ,x) - SS(𝒫,x))/𝒫.σ2
	Υ = Δll  + log(αᵒ) - log(α) + logpdf(prior_α, αᵒ) - logpdf(prior_α, α)
	if log(rand()) < Υ
		α = αᵒ
		acc[1] += 1
	end
	α, acc
end

"""
	update_ξ(𝒫,x)
"""
function update_ξ(𝒫, x, σ2γ)
	n = length(x)
	ᾱ = [sqrt(𝒫.Δ[i]) * 𝒫.α / (1.0 + 𝒫.α* 𝒫.Δ[i]) for i ∈ eachindex(𝒫.Δ)]
	U = ec1([(x[i] - A(i-1,𝒫) * x[i-1])/sqrt(𝒫.Δ[i-1]) for i ∈ 2:n])
	nb = nbasis(𝒫.J)
	V = zeros(nb, nb)
	v = zeros(nb)
	for i ∈ 2:n
		#global v, V
		ϕi = ϕ(𝒫.t[i], 𝒫.J)
		V += ᾱ[i-1]^2 * ϕi * ϕi'
		v += ᾱ[i-1] * U[i-1] * ϕi
	end
	V = PDMat(Symmetric(V/𝒫.σ2 + I/σ2γ))
	return rand(MvNormalCanon(v/𝒫.σ2,V))
end


function mcmc(𝒫, y; ITER = 1000,
				propσ_α=0.2, propσ_ψ=0.2,
	 			prior_α = Exponential(10.0),
				prior_ψ1 = Uniform(0,1),
				Aσ=0.1, Bσ=0.1,
				Aψ=0.1, Bψ=0.1,
				σ2γ = 10.0,
				print_skip=100)

	#m0= zeros(1) ;

	m0 = [mean(ec1(y))]
	P0 = 0.1*Matrix(1.0I, d, d)
	𝒢 = ObsGroup(𝒫, y)
	θ = [parameters(𝒫)]
	X = []
	acc = [0,0]
	for it ∈ 2:ITER

		if mod(it,print_skip)==0 println("iteration $it") end
		(m, P), (m⁻, P⁻) = ff(y, (m0,P0), 𝒫)
		xs = bsample(y, (m, P), (m⁻, P⁻), 𝒫)

		ψ, acc  = update_ψ(𝒫, 𝒢, xs, acc, prior_ψ1, Aψ, Bψ, propσ_ψ)
		σ2 = update_σ2(𝒫, xs, Aσ, Bσ)
		𝒫 = DF(𝒫.α, 𝒫.ξ,  σ2, ψ, 𝒫.t, 𝒫.Δ, 𝒫.typeobs, 𝒫.J, 𝒫.η)

		α, acc = update_α(𝒫, xs, acc, prior_α, propσ_α)
		𝒫 = DF(α, 𝒫.ξ, 𝒫.σ2, 𝒫.ψ, 𝒫.t, 𝒫.Δ, 𝒫.typeobs, 𝒫.J, 𝒫.η)

		ξ= update_ξ(𝒫, xs, σ2γ)
		𝒫 = DF(𝒫.α, ξ, 𝒫.σ2, 𝒫.ψ, 𝒫.t, 𝒫.Δ, 𝒫.typeobs, 𝒫.J, 𝒫.η)

		push!(θ, parameters(𝒫))
		push!(X, deepcopy(xs))
	end
	accperc_α = round.(100acc/(ITER-1);digits=2)
	println("Acceptance percentage for updating α: $accperc_α%")
	θ, X, 𝒫, accperc_α
end


ec(x,i) = map(u->u[i],x)
ec1(x) = ec(x,1)
