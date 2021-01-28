"""
	readdata(path)

Processing of input data.

Returns:
t:: vector containing observation times
Δ:: vector with lengths of interobservation time intervals
y:: vector of static vectors, containing the measured chorophyl-A concentrations
typeobs:: vector of characters containing info on type of observation at each time instance
(obs1: only in situ measurment; obs2: only satellite measurment; obs3: both measurements)
"""
function readdata(path)
	dat = CSV.File(path) |> DataFrame!
	t = vcat(0.0,dat[!,:time_elapsed])
	typeobs = String.(dat.obsscheme)#dat[!,:obsscheme]
	Δ = diff(t)
	y = SVector[]
	for r in eachrow(dat)
	    if r[:obsscheme]=="obs1"
	        u = r[:chl_water]
	        push!(y,SVector(u))
	    elseif r[:obsscheme]=="obs2"
	        u = r[:chl]
	        push!(y,SVector(u))
	    elseif r[:obsscheme]=="obs3"
	        u = [r[:chl_water], r[:chl]]
	        push!(y,SVector{2}(u))
	    end
	end
	lrad_temp = dat[!,[:lograd, :temp]]
	t, Δ, y, typeobs, Matrix(lrad_temp)
end


# Model: d X_t = - α (X_t - μ(t)) d t + σ d W_t
# μ(t) = ∑ᵢ ξᵢ
# extrinsic noise with variances governed by ψ

struct DF{T}
    α::Float64  # mean reversion par
    ξ::Vector{Float64}  #  pars in periodic drift function
    σ2::Float64    # squared diffusivity
	ψ̄::Float64
    ψ::Float64  # vars on observation equation
    t::Vector{Float64}  # t_i i ∈ 1...n
    Δ::Vector{Float64}  # Δ_i = t_i-t_{i-1}
    typeobs::Vector{T}
	J::Int64			# truncation level Fourier series (so 2J+1 basis functions)
	η::Float64
	lrad_temp::Matrix{Float64}  # first column contains lograd; second column temperature
	ρ::Vector{Float64}  # pars for including logradiation and temperature
end

parameters(𝒫::DF) = (𝒫.α, 𝒫.ξ, 𝒫.σ2, 𝒫.ψ̄, 𝒫.ψ, 𝒫.ρ)

struct ObsGroup{T}
	ind1::Vector{Int64} # indices in y where measurement device 1 is used
	ind2::Vector{Int64} # indices in y where measurement device 2 is used
	y1::Vector{T} 			# y where measurement device 1 is used
	y2::Vector{T}			# y where measurement device 2 is used

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
		new{eltype(y1)}(ind1,ind2,y1,y2)
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
nbasis(J) = 2^(J)



μ(t,ξ,J) = dot(ξ, ϕ(t, J))


if implicit
    A(k,𝒫) = k==0 ? (@SMatrix [1.0]) : @SMatrix [(1.0 + 𝒫.α * 𝒫.Δ[k])^(-1)]
    a(k,𝒫) = k==0 ? SVector(0.0) : SVector((1.0 + 𝒫.α * 𝒫.Δ[k])^(-1) * 𝒫.α * μ(𝒫.t[k+1],𝒫.ξ, 𝒫.J) * 𝒫.Δ[k])
else
    A(k,𝒫) = k==0 ?  1.0 : 1.0 - 𝒫.α * 𝒫.Δ[k]
    a(k,𝒫) = k==0 ?  0.0 : 𝒫.α *  μ(𝒫.t[k+1],𝒫.ξ) * 𝒫.Δ[k]
end

Q(k,𝒫) = k==0 ?  (@SMatrix [0.0]) :  @SMatrix [𝒫.σ2 * 𝒫.Δ[k]]

H(k,𝒫) =  𝒫.typeobs[k]=="obs3" ? (@SMatrix [1.0; 1.0]) : (@SMatrix [1.0])

function R(k,𝒫)
	if   𝒫.typeobs[k]=="obs1"
    	return @SMatrix [𝒫.η * 𝒫.ψ̄ * 𝒫.ψ]
    elseif    𝒫.typeobs[k]=="obs2"
        return @SMatrix  [𝒫.ψ]
    else
        return @SMatrix  [𝒫.η * 𝒫.ψ̄ * 𝒫.ψ 0.0; 0.0 𝒫.ψ]
    end
end


## mcmc updates for pars

function SS(𝒫, x)
	S = 0.0
	for k ∈ 2:length(x)
		#S += (x[k] - A(k-1,𝒫) * x[k-1] - a(k-1,𝒫)).^2/𝒫.Δ[k-1]
		S += norm(x[k] - A(k-1,𝒫) * x[k-1] - a(k-1,𝒫))^2/𝒫.Δ[k-1]
	end
	S
end

function update_σ2(𝒫, x, Aσ, Bσ)
	m = length(x)-1
	S = SS(𝒫, x)
	rand(InverseGamma(0.5m+ Aσ, 0.5S + Bσ))
end

logtargetψ̄(ψ̄, ψ ,n ,S, η) = -0.5n * log(ψ̄) - 0.5S/(η * ψ̄ * ψ)  # needed for MH-update first element of ψ

function update_ψ̄ψ(𝒫::DF, 𝒢::ObsGroup, x, acc, prior_ψ̄, Aψ, Bψ, propσ_ψ̄)
	S1 = norm(ec1(𝒢.y1-x[𝒢.ind1]))^2
	S2 = norm(ec1(𝒢.y2-x[𝒢.ind2]))^2
	n1 = length(𝒢.y1)
	n2 = length(𝒢.y2)
	ψ = rand(InverseGamma(Aψ + 0.5*(n1 + n2), Bψ + 0.5*(S1/(𝒫.η * 𝒫.ψ̄) + S2) ))
	ψ̄ = 𝒫.ψ̄
	ψ̄ᵒ = ψ̄/(ψ̄ + (1-ψ̄)* exp(propσ_ψ̄ * randn()))
	Δll = logtargetψ̄(ψ̄ᵒ,ψ,n1,S1,𝒫.η) - logtargetψ̄(ψ̄,ψ,n1,S1,𝒫.η)
	log_jacob_term = log(ψ̄ᵒ*(1.0-ψ̄ᵒ)) - log(ψ̄*(1.0-ψ̄))
	Υ = Δll + log_jacob_term  + logpdf(prior_ψ̄, ψ̄ᵒ) - logpdf(prior_ψ̄, ψ̄)
	if log(rand()) < Υ
		ψ̄ = ψ̄ᵒ
		acc[2] += 1
	end
	return ψ̄, ψ, acc
end


function update_α(𝒫, x, acc, prior_α, propσ)
	α = 𝒫.α
	αᵒ = α * exp(propσ*randn())
	𝒫ᵒ = DF(αᵒ, 𝒫.ξ, 𝒫.σ2, 𝒫.ψ̄, 𝒫.ψ, 𝒫.t, 𝒫.Δ, 𝒫.typeobs, 𝒫.J, 𝒫.η, 𝒫.lrad_temp, 𝒫.ρ)
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
function update_ξ(𝒫, x, priorvarξρ)
	n = length(x)
	ᾱ = [sqrt(𝒫.Δ[i]) * 𝒫.α / (1.0 + 𝒫.α* 𝒫.Δ[i]) for i ∈ eachindex(𝒫.Δ)]
	U = ec1([(x[i] - A(i-1,𝒫) * x[i-1])/sqrt(𝒫.Δ[i-1]) for i ∈ 2:n])
	nb = nbasis(𝒫.J)
	Z = zeros(nb, nb)
	v = zeros(nb)
	for i ∈ 2:n
		zi = ϕ(𝒫.t[i], 𝒫.J) * ᾱ[i-1]
		Z += zi * zi'
		v += U[i-1] * zi
	end
	P = PDMat(Symmetric(Z/𝒫.σ2 + I/priorvarξρ))
	ν = v/𝒫.σ2
	rand(MvNormalCanon(ν,P))
end

"""
	mcmc

printskip: output every print_skip iteration
saveskip: nunber of sampled paths to skip in saving to output
"""
function mcmc(𝒫, y; ITER = 1000,
				propσ_α=0.2, propσ_ψ̄=0.2,
	 			prior_α = Exponential(10.0),
				prior_ψ̄ = Uniform(0,1),
				Aσ=0.1, Bσ=0.1,
				Aψ=0.1, Bψ=0.1,
				priorvarξ = 10.0, # prior var on ξ_i and ρ_j
				printskip=50,
				saveskip=25
				)
	m0 = @SVector [0.0]
	P0 = @SMatrix [10.0]

	𝒢 = ObsGroup(𝒫, y)
	θ = [parameters(𝒫)]
	acc = [0,0]

	(m, P), (m⁻, P⁻) = ff(y, (m0,P0), 𝒫)
	xs = bsample((m, P), (m⁻, P⁻), 𝒫)
	X = [xs]

	for it ∈ 2:ITER

		if it % printskip == 0 println("iteration $it") end

		ψ̄, ψ, acc  = update_ψ̄ψ(𝒫, 𝒢, xs, acc, prior_ψ̄, Aψ, Bψ, propσ_ψ̄)
		𝒫 = @set 𝒫.ψ̄ = ψ̄
		𝒫 = @set 𝒫.ψ = ψ

		σ2 = update_σ2(𝒫, xs, Aσ, Bσ)
		𝒫 = @set 𝒫.σ2 = σ2

		α, acc = update_α(𝒫, xs, acc, prior_α, propσ_α)
		𝒫 = @set 𝒫.α = α

		ξ = update_ξ(𝒫, xs, priorvarξ)
		𝒫 = @set 𝒫.ξ = ξ

		ff!(y, (m0,P0), (m, P), (m⁻, P⁻), 𝒫)
		bsample!(xs, (m, P), (m⁻, P⁻), 𝒫)

		push!(θ, parameters(𝒫))
		if it % saveskip == 0
			push!(X, deepcopy(xs))
		end
	end
	accperc_α = round.(100acc/(ITER-1);digits=2)
	println("Acceptance percentage for updating (α, ψ̄): $accperc_α%")

	θ, X, 𝒫, accperc_α
end


ec(x,i) = map(u->u[i],x)
ec1(x) = ec(x,1)
