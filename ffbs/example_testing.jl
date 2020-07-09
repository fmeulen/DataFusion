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



α = 10.0
J = 5
ξ = vcat([0.2, 1.0, 0.0], zeros(nbasis(J)-3))
σ2 = .1^2
ψ = [.01, .01]

N = 150
Δ = rand(Uniform(0.1,0.15),N)
t = vcat(0.0,cumsum(Δ))
typeobs = vcat(fill("obs1",N-20), fill("obs3",10),fill("obs2",10))
𝒫true = DF(α,  ξ,  σ2, ψ, t, Δ, typeobs,J)

d = size(A(1,𝒫true))[1]
x0 = zeros(d)
x = [x0]
y =[] # irrelevant but of right type # observations are at indices 2:N
for k in 1:N
    global x
    xk = k==1 ? x0 : A(k-1,𝒫true)*x[k-1] + a(k-1,𝒫true) + rand(Gaussian(zeros(d), Q(k-1,𝒫true)))
    push!(x, xk)
    yk =  rand(Gaussian( H(k,𝒫true)*xk[1], R(k,𝒫true)))
    push!(y,SVector(yk))
end

pl = Plots.plot(t[2:end],ec1(y))
m0= zeros(d) ; P0=0.0*Matrix(1.0I, d, d)
(m, P), (m⁻, P⁻) = ff(y, (m0,P0), 𝒫true)
M = 3 # nr of ffbs pats
xs = [bsample(y, (m, P), (m⁻, P⁻), 𝒫true) for _ ∈ 1:M]
for k ∈ 1:M
    Plots.plot!(pl,t[2:end],ec1(xs[k]))
end
display(pl)
Plots.plot!(pl,t[2:end],ec1(m))

####################################################

# initialise 𝒫
α = 5.8
ξ = vcat([0.5, 0.0, 1.0], rand(nbasis(J)-3))
σ2 = 0.8
ψ = [0.08, 0.08]
𝒫init = DF(α,  ξ,  σ2, ψ, t, Δ, typeobs, J)

ITER = 5000
θ, X, 𝒫, accperc_α = mcmc(𝒫init, y; ITER = ITER , propσ=0.1)

include("postprocessing.jl")
