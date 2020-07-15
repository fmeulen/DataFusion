using DataFrames
using Kalman, GaussianDistributions, LinearAlgebra
using GaussianDistributions: ⊕
using Plots
using Distributions
using StaticArrays
using CSV
using RCall
using PDMats

cd(@__DIR__)
include("ffbs.jl")
include("funcdefs.jl")

## Lorinc example: read true data
α = 14.0
J = 5
ξ = vcat([0.5, 1.0, 2.0], zeros(nbasis(J)-3))
σ2 = .1^2
ψ = [.01, .01] # first el should be <1
η = 658/8005


dat_all = CSV.read("observations.csv")
K = 8005;   dat = dat_all[1:K,:]
t = vcat(0.0,dat[!,:time_elapsed])
typeobs = dat[!,:obsscheme]
Δ = diff(t)
𝒫 = DF(α, ξ, σ2, ψ, t, Δ, typeobs, J, η)

trf(x) = log(x) # apply model to the log of the data
invtrf(x) = exp(x)

y = []
for r in eachrow(dat)
	if r[:obsscheme]=="obs1"
		u = tryparse(Float64,r[:chl_water])
		push!(y,SVector(trf(u)))
	elseif r[:obsscheme]=="obs2"
		u = tryparse(Float64,r[:chl])
		push!(y,SVector(trf(u)))
	elseif r[:obsscheme]=="obs3"
		u = tryparse.(Float64,[r[:chl_water], r[:chl]])
		push!(y,SVector{2}(trf.(u)))
	end
end

#row 7627 is of type obs3  r= dat[7627,:]

d = 1
m0= zeros(d) ; P0=0.0*Matrix(1.0I, d, d)
(m, P), (m⁻, P⁻) = ff(y, (m0,P0), 𝒫)
pl = Plots.plot(t[2:end],first.(y))
M = 5 # nr of ffbs pats
Ys = [bsample(y, (m, P), (m⁻, P⁻), 𝒫) for _ ∈ 1:M]
for k ∈ 1:M
    Plots.plot!(pl,t[2:end],ec1(Ys[k]))
end
display(pl)



####################################################

# initialise 𝒫
𝒫init = DF(α,  ξ,  σ2, ψ, t, Δ, typeobs, J, η)

# apriori expect water measurement to be more accurate, i.e. ψ1 smaller than ψ2

ITER = 1000
θ, X, 𝒫, accperc_α = mcmc(𝒫init, y; ITER = ITER , propσ_α=0.1, propσ_ψ=0.2)

𝒫true = 𝒫init # simply unknown here
include("postprocessing.jl")
