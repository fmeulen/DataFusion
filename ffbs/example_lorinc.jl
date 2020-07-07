using DataFrames
using Kalman, GaussianDistributions, LinearAlgebra
using GaussianDistributions: ⊕
using Plots
using Distributions
using StaticArrays
using CSV
using RCall

cd(@__DIR__)
include("ffbs.jl")
include("funcdefs.jl")

## Lorinc example: read true data

dat_all = CSV.read("observations.csv")
K = 8005;   dat = dat_all[1:K,:]
t = vcat(0.0,dat[:time_elapsed])
typeobs = dat[:obsscheme]
Δ = diff(t)
𝒫 = DF(α, ξ, c, σ2, ψ, t, Δ, typeobs)

y = []
for r in eachrow(dat)
	if r[:obsscheme]=="obs1"
		u = tryparse(Float64,r[:chl_water])
		push!(y,SVector(u))
	elseif r[:obsscheme]=="obs2"
		u = tryparse(Float64,r[:chl])
		push!(y,SVector(u))
	elseif r[:obsscheme]=="obs3"
		u = tryparse.(Float64,[r[:chl_water], r[:chl]])
		push!(y,SVector{2}(u))
	end
end

#row 7627 is of type obs3
r= dat[7627,:]

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
α = .8
c = 5.0#𝒫true.c#5.0
ξ = [6.0, 2.0]
σ2 = 0.8
ψ = [0.08, 0.08]
𝒫init = DF(α,  𝒫true.ξ,  c, σ2, ψ, t, Δ, typeobs)

ITER = 1000
θ, X, 𝒫, accperc_α = mcmc(𝒫init, y; ITER = ITER , propσ=0.2)

𝒫true = 𝒫init # simply unknown here
include("postprocessing.jl")
