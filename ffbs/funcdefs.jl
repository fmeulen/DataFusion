struct DF
    α::Float64  # mean reversion par
    ξ::Vector{Float64}  #  pars in periodic drift function
	c::Float64     # multiplicative constant in drift function
    σ2::Float64    # squared diffusivity
    ψ::Vector{Float64}  # vars on observation equation
    t::Vector{Float64}  # t_i i ∈ 1...n
    Δ::Vector{Float64}  # Δ_i = t_i-t_{i-1}
    typeobs::Vector
end

parameters(𝒫::DF) = (𝒫.α, 𝒫.ξ, 𝒫.c, 𝒫.σ2, 𝒫.ψ)

struct ObsGroup
	ind1::Vector{Int64} # indices in y where measurement device 1 is used
	ind2::Vector{Int64} # indices in y where measurement device 2 is used
	y1::Vector 			# y where measurement device 1 is used
	y2::Vector			# y where measurement device 2 is used
end

# Model equations
# x[k] = A(k-1) * x[k-1] + a[k-1] +  N(0,Q[k-1])
# y[k] = H(k) * x[k] + N(0,R(k))

# We assume typeobs="obs1" (only meas device 1 used), typeobs="obs2" (only meas device 2 used),
# typeobs="obs3" (both meas devices 1 and 2 used)

implicit = true# true
μ(t,ξ) =          pdf(Beta(ξ[1],ξ[2]),mod(t,1.0))#dot(ξ, ϕ(t))
if implicit
    A(k,𝒫) = k==0 ?  SMatrix{1,1}([1.0]) : SMatrix{1,1}( [(1.0 + 𝒫.α * 𝒫.Δ[k])^(-1)] )
    a(k,𝒫) = k==0 ?  (@SVector [0.0]) :   (@SVector [A(k,𝒫)[1,1] * 𝒫.α * 𝒫.c * μ(𝒫.t[k+1],𝒫.ξ) * 𝒫.Δ[k]   ])
else
    A(k,𝒫) = k==0 ?  SMatrix{1,1}([1.0]) : SMatrix{1,1}( [1.0 - 𝒫.α * 𝒫.Δ[k]] )
    a(k,𝒫) = k==0 ?  (@SVector [0.0]) :   (@SVector [𝒫.α * 𝒫.c * μ(𝒫.t[k+1],𝒫.ξ) * 𝒫.Δ[k]   ])
end
Q(k,𝒫) = k==0 ?  SMatrix{1,1}([0.0]) :   SMatrix{1,1}( [𝒫.σ2 * 𝒫.Δ[k]] )
H(k,𝒫) =  𝒫.typeobs[k]=="obs3" ? SMatrix{2,1}([1.0 1.0]) : SMatrix{1,1}([1.0])
function R(k,𝒫)
    if   𝒫.typeobs[k]=="obs1"
        return  SMatrix{1,1}([𝒫.ψ[1]])
    elseif    𝒫.typeobs[k]=="obs2"
        return  SMatrix{1,1}([𝒫.ψ[2]])
    else
        return  SMatrix{2,2}([𝒫.ψ[1] 0.0; 0.0 𝒫.ψ[2]])
    end
end

function grouping(𝒫, y)
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
	ObsGroup(ind1,ind2,y1,y2)
end

## mcmc updates for pars
function SS(𝒫, x)
	S = 0.0
	for k ∈ 2:length(x)
		S += (x[k] - A(k-1,𝒫) * x[k-1] - a(k-1,𝒫))[1]^2/𝒫.Δ[k-1]
	end
	S
end

function update_σ2(𝒫, x; Aσ=.01, Bσ=0.01)
	m = length(x)-1
	S = SS(𝒫, x)
	rand(InverseGamma(0.5m+ Aσ, 0.5S + Bσ))
end

function update_ψ(𝒫::DF, 𝒢::ObsGroup, x; Aσ=0.01, Bσ=0.01)
	S1 = norm(ec1(𝒢.y1-x[𝒢.ind1]))^2
	S2 = norm(ec1(𝒢.y2-x[𝒢.ind2]))^2
	ψ1 = rand(InverseGamma(Aσ + 0.5length(𝒢.y1), Bσ + 0.5S1 ))
	ψ2 = rand(InverseGamma(Aσ + 0.5length(𝒢.y2), Bσ + 0.5S2 ))
	return [ψ1, ψ2]
end

function update_c(𝒫,x)
	𝒫1 = DF(𝒫.α, 𝒫.ξ, 1.0, 𝒫.σ2, 𝒫.ψ, 𝒫.t, 𝒫.Δ, 𝒫.typeobs)  # 𝒫 with c=1
	S1 = 0.0
	S2 = 0.0
	for k in 2:length(x)
		S1 += ((x[k] - A(k-1,𝒫1) * x[k-1]) .* a(k-1,𝒫1)./Q(k-1,𝒫1))[1,1]
		S2 += a(k-1,𝒫1)[1,1]^2/Q(k-1,𝒫1)[1,1]
	end
	rand(Normal(S1/S2, 1/√S2))
end

function update_α(𝒫, x, acc ; propσ=0.1)
	α = 𝒫.α
	αᵒ = α * exp(propσ*randn())
	𝒫ᵒ = DF(αᵒ, 𝒫.ξ, 𝒫.c, 𝒫.σ2, 𝒫.ψ, 𝒫.t, 𝒫.Δ, 𝒫.typeobs)
	Δll = -0.5*(SS(𝒫ᵒ,x) - SS(𝒫,x))/𝒫.σ2
	if log(rand()) <  (Δll  + log(αᵒ) - log(α)) #+ logpdf(prior[1],λᵒ) - logpdf(prior[1],P.λ)
		α = αᵒ
		acc += 1
	end
	α, acc
end

function mcmc(𝒫, y; ITER = 1000, propσ=0.2)
	m0= zeros(d) ; P0=0.0*Matrix(1.0I, d, d)
	𝒢 = grouping(𝒫, y)
	θ = [parameters(𝒫)]
	X = []
	acc = 0
	for it ∈ 2:ITER

		if mod(it,50)==0 println("iteration $it") end
		(m, P), (m⁻, P⁻) = ff(y, (m0,P0), 𝒫)
		xs = bsample(y, (m, P), (m⁻, P⁻), 𝒫)

		ψ = update_ψ(𝒫, 𝒢, xs)
		σ2 = update_σ2(𝒫, xs)
		𝒫 = DF(𝒫.α, 𝒫.ξ, 𝒫.c, σ2, ψ, 𝒫.t, 𝒫.Δ, 𝒫.typeobs)

		α, acc = update_α(𝒫, xs, acc; propσ=propσ)
		𝒫 = DF(α, 𝒫.ξ, 𝒫.c, 𝒫.σ2, 𝒫.ψ, 𝒫.t, 𝒫.Δ, 𝒫.typeobs)

		c = update_c(𝒫, xs)#𝒫true.c
		𝒫 = DF(𝒫.α, 𝒫.ξ, c, 𝒫.σ2, 𝒫.ψ, 𝒫.t, 𝒫.Δ, 𝒫.typeobs)

		push!(θ, parameters(𝒫))
		push!(X, deepcopy(xs))
	end
	accperc_α = round(100acc/(ITER-1);digits=2)
	println("Acceptance percentage for updating α: $accperc_α%")
	θ, X, 𝒫, accperc_α
end


ec(x,i) = map(u->u[i],x)
ec1(x) = map(u->u[1],x)



# function ν(x,a0,a,b)
# 	N = length(a)
# 	S = a0
# 	for  n ∈ eachindex(b)
# 		S += a[n]cos(2π*n*x) + b[n]sin(2π*n*x)
# 	end
# 	return S
# end
# ν(a0,a,b) = x -> ν(x,a0,a,b)
#
# x = range(-2.0, 1.0; length=100)
# plot(x,ν(1.0,randn(6),randn(6)).(x))



function ϕ(x; J=5)
	out = [1.0]
	for j ∈ 1:J
		push!(out, cos(2π*j*x))
		push!(out, sin(2π*j*x))
	end
	out
end

# function ϕ(x,k)
# 	if !isinteger(k) error("k should be integer valued.") end
# 	if k==1 return 1.0
# 	elseif isodd(k) return sin(π*x*(k-1)) #sin(2π*x*0.5(k-1))
# 	else return cos(π*x*k) #cos(2π*x*0.5k)
# 	end
# end
# ϕ(k) = x-> ϕ(x,k)

"""
	update_ξ(𝒫,x)
"""
function update_ξ(𝒫, x)
	n = length(x)
	ᾱ = [𝒫.α/(1.0 + 𝒫.α* 𝒫.Δ[i]) for i ∈ eachindex(𝒫.Δ)]
	U = ec1([(x[i] - A(i-1,𝒫) * x[i-1])/𝒫.Δ[i-1] for i ∈ 2:n])
	V = zeros(2J+1, 2J+1)
	v = zeros(2J+1)
	for i ∈ 2:n
		ϕi = ϕ(x[i-1][1,1])
		V += ᾱ[i-1]^2 * ϕi * ϕi'
		v += ᾱ[i-1] * U[i-1] * ϕi
	end
	V = PDMat(Symmetric(V))
	rand(MvNormalCanon(v,V))
end
