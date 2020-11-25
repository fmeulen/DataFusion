using DataFrames
using GaussianDistributions, LinearAlgebra
#using Plots
using Distributions
using StaticArrays
using CSV
#using RCall
using PDMats
using DelimitedFiles

cd(@__DIR__)
include("ffbs.jl")
include("fdefs.jl")

## Lorinc example: read true data, ensure log-chlorA, log-radiation and temperature are given
t, Δ, conc, typeobs, lrad_temp = readdata("../../csv/observation_lograd_temp.csv")
# change to logconcentration
logconc = map( x-> log.(x), conc)
####################################################
α = 14.0
J = 5
ξ = vcat([0.5, 1.0, 2.0], zeros(nbasis(J)-3))
σ2 = .1^2
ψ̄ = 0.01 # should be <1
ψ = 0.01
η = 658/8005
ρ = [0.0, 0.0]

####################################################

# initialise 𝒫
𝒫init = DF(α,  ξ,  σ2, ψ̄, ψ, t, Δ, typeobs, J, η, lrad_temp, ρ)

# apriori expect water measurement to be more accurate, i.e. ψ1 smaller than ψ2

ITER = 100
θ, X, 𝒫, accperc_α, postmean_paths = mcmc(𝒫init, logconc; ITER = ITER)

## postprocessing

writedlm("../../csv/postmean_paths.csv", postmean_paths)

using Plots
### traceplots
# parameters(𝒫::DF) = (𝒫.α, 𝒫.ξ, 𝒫.σ2, 𝒫.ψ̄, 𝒫.ψ, 𝒫.ρ)
p1 = plot( ec(map(x->x[1],θ),1),label="α")
p2 = plot( ec(map(x->x[2],θ),1), label="ψ1")
p3 = plot( ec(map(x->x[3],θ),1), label="σ²")
p4 = plot( ec(map(x->x[4],θ),1), label="ψ̄")
p5 = plot( ec(map(x->x[5],θ),1), label="ψ")
p6 = plot( ec(map(x->x[6],θ),1), label="ρ₁")
p6 = plot( ec(map(x->x[6],θ),2), label="ρ₂")

l = @layout [a b c ; d e f]
plot(p1, p2, p3, p4, p5, p6, layout = l)

pX1 = plot(ec1(X[1]))
pXhalf = plot(ec1(X[div(ITER,2)]))
pXend = plot(ec1(X[ITER-1]))
plot(pX1, pXhalf, pXend, layout=@layout [a;b;c])


BI = div(ITER,2)
using RCall
include("postprocessing.jl")
