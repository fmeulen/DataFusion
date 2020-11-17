using DataFrames
using GaussianDistributions, LinearAlgebra
using Plots
using Distributions
using StaticArrays
using CSV
#using RCall
using PDMats
using DelimitedFiles

cd(@__DIR__)
include("ffbs.jl")
include("funcdefs.jl")

## Lorinc example: read true data, ensure log-chlorA, log-radiation and temperature are given
t, Δ, conc, typeobs, lrad_temp = readdata("../../csv/observation_lograd_temp.csv")
# change to logconcentration
logconc = map( x-> log.(x), conc)
####################################################

# initialise 𝒫
α = 14.0
J = 5
ξ = vcat([0.5, 1.0, 2.0], zeros(nbasis(J)-3))
σ2 = .1^2
ψ̄ = 0.01
ψ = 0.01


η = [0.0, 0.0]
𝒫init = DF(α,  ξ,  σ2, ψ̄, ψ, t, Δ, typeobs, J, η, lrad_temp)

# apriori expect water measurement to be more accurate, i.e. ψ1 smaller than ψ2
ITER = 1000
θ, X, 𝒫, accperc_α, postmean_paths = mcmc(𝒫init, conc; ITER = ITER, propσ=0.1)

writedlm("../../csv/postmean_paths.csv", postmean_paths)

using Plots
p = plot(ec1(X[1]))
plot!(p, ec1(X[ITER-1]))

BI = div(ITER,2)
using RCall
include("postprocessing.jl")

# plot η1, η2
pp = plot( ec(map(x->x[5],θ),1))
plot!(pp, ec(map(x->x[5],θ),2))
