using Revise
using GuidedProposals
using DiffusionDefinition, ObservationSchemes
const GP = GuidedProposals
const DD = DiffusionDefinition
const OBS = ObservationSchemes

using StaticArrays, Random, Plots
using LinearAlgebra
using Distributions
using CSV
# using PyPlot
# pyplot()

# ᾱ = log α, σ̄ = log σ

workdir = @__DIR__
cd(workdir)
include("funcdefs.jl")
include("model.jl")

#---------------------------------------------------------
# ordering of parameters is (ᾱ, ω̄, ϕ̄, δ, σ̄)
#---------------------------------------------------------

# Data generation
θtrue = [log(6.0), log(2.0), log(5.0), 0.8, log(0.15)]
Ptrue = Periodic1D(θtrue...)
tt, y1 = 0.0:0.001:8.0, @SVector [1.0]
X = rand(Ptrue, tt, y1)
data = map(
	x->(x[1], x[2][1] + 0.1randn()),
	collect(zip(X.t, X.x))[1:75:end]
)[2:end]
# let's examine the data
gr()
Plots.plot(X, Val(:vs_time))#, size=(800, 300))
scatter!(map(x->x[1], data), map(x->x[2], data), label="data")

θinit = Dict(:ᾱ=>log(1.0), :ω̄=>log(3.3), :ϕ̄=> log(4.0), :δ=>1.0, :σ̄=>log(3.0))
θinit_full = values(θinit)

# extrinsic noise
σ1 = 3.0
σ2 = 2.0 # should be bigger than σ1

K = 8005 # 270
include("setobs_adv.jl")

# plot discrete time observations
p = Plots.plot()#size=(1400, 800))
Plots.scatter!(p,ts[1:K], map(x->x[1], X[1:K]),color="red",label="data",markersize=2)
display(	p)
#png(p,"data.png")

# p = Plots.plot()#size=(1400, 800))
# Plots.scatter!(p,mod.(ts[1:K],1), map(x->x[1], X[1:K]),color="red",label="data",markersize=2)
# display(	p)

τ(t₀,T) = (t) ->  t₀ + (t-t₀) * (2-(t-t₀)/(T-t₀))
standard_guid_prop_time_transf(tt) = τ(tt[1], tt[end]).(tt)


# rec = recording.obs[1]
dt = 0.1
grids = OBS.setup_time_grids(recording, dt)

impgrid = Array{Float64,1}[]
nimp = vcat(fill(15, 270), fill(5,8005-270))
tt1 = collect(range(recording.t0,recording.obs[1].t,length=nimp[1]))
push!(impgrid, standard_guid_prop_time_transf(tt1))
for k ∈ 2:8005
	tt1 = collect(range(recording.obs[k-1].t,recording.obs[k].t,length=nimp[k]))
	push!(impgrid, standard_guid_prop_time_transf(tt1))
end
length(recording.obs)
#length(impgrid)-length(grids)
impgrid

#PP = build_guid_prop(Periodic1DAux, recording, grids)
PP = build_guid_prop(Periodic1DAux, recording, impgrid)

@time 	paths, θθ, Q = simple_inference(
	Periodic1DAux, recording, impgrid, θinit;
	sc =0.05, ρ=0.01, num_steps=200, λ_cov = 0.1, λ_Id = 0.1)

p = Plots.plot(size=(1400, 800))
pathsplot = paths[1:10]
pathsplot = paths[end-2:end]
for path in pathsplot,  i in eachindex(path)
		Plots.plot!(p, path[i], Val(:vs_time), alpha=0.4, label="", color=["steelblue"])
end
Plots.scatter!(p,ts[1:K], map(x->x[1], X[1:K]),color="red",label="data",markersize=8)
display(p)
png(p,"reconstruction.png")
println(θθ[3])


# Examine the results
println(DD.parameters(Ptrue))

# ordering ᾱ, ω1, σ̄, δ1 in the following plot
#Plots.plot(θθ, labels=["ᾱ", "ω̄", "ϕ̄", "δ", "σ̄"])

p1 = Plots.plot(θθ[1],label="ᾱ")
#hline!(p1, [DD.parameters(P)[:ᾱ]],label="")
p2 = Plots.plot(θθ[2],label="ω̄")
#hline!(p2, [DD.parameters(P)[:ω̄]],label="")
p3 = Plots.plot(θθ[3], label="ϕ̄")
#hline!(p3, [DD.parameters(P)[:ϕ̄]],label="")
p4 = Plots.plot(θθ[4], label="δ")
#hline!(p4, [DD.parameters(P)[:δ]],label="")
p5 = Plots.plot(θθ[5], label="σ̄")
#hline!(p5, [DD.parameters(P)[:σ̄]],label="")

Plots.plot(p1,p2,p3,p4, p5)


if false
	p0 = Plots.plot(size=(1400, 800))
	scatter!(p0,map(x->x[1], data), map(x->x[2], data), label="data", markersize=8, markercolor=:grey)
	display(p0)
	png(p0,"data.png")

	p = Plots.plot(size=(1400, 800))
	for path in paths[end-9:end]
		for i in eachindex(path)
			Plots.plot!(p, path[i], Val(:vs_time), alpha=0.4, label="", color=["steelblue"])
		end
	end
	Plots.scatter!(p,ts[1:K], map(x->x[1], X[1:K]),color="red",label="data")

	# Plots.plot!(p,X, Val(:vs_time),color="red", alpha=0.4)
	# scatter!(p, map(x->x[1], data), map(x->x[2], data), label="data", markersize=8, markercolor=:grey)
	display(p)
	png(p,"reconstruction.png")
end
