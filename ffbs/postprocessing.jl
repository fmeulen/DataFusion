using DelimitedFiles

BI = div(ITER,2)
postmean_paths = ec1([mean(map(x->x[i],X[BI:ITER-1])) for i in eachindex(X[1])])


pl = Plots.scatter(𝒫.t[2:end],first.(y), markersize=3)
for k ∈ eachindex(X[end-10:end])
    Plots.plot!(pl,𝒫.t[2:end],ec1(X[k]),label="")
end
Plots.plot!(pl, 𝒫.t[2:end],postmean_paths)
display(pl)
#png(pl,"test1.png")
writedlm("/Users/Frank/.julia/dev/DataFusion/processing_in_r/postmean_paths.csv",postmean_paths)


p1=Plots.plot( ec(θ,1),label="α")
p2=Plots.plot( first.(ec(θ,2)),label="ξ1")
p3=Plots.plot( ec(θ,3),label="σ2")
p4=Plots.plot( first.(ec(θ,4)),label="ψ1")
p5=Plots.plot( last.(ec(θ,4)),label="ψ2")
#png(p, "pars.png")


df = DataFrame(iterate= repeat(1:ITER,6),
	parameter= vcat(ec(θ,1),ec(θ,3),first.(ec(θ,4)),last.(ec(θ,4)),first.(ec(θ,2)),last.(ec(θ,2))),
	type=repeat(["alpha","sigma2","psi1","psi2","xi1","xilast"],inner=ITER))
dftrue = DataFrame(type=["alpha","sigma2","psi1","psi2","xi1","xilast"], parameter=[𝒫true.α, 𝒫true.σ2, 𝒫true.ψ[1], 𝒫true.ψ[2], 𝒫true.ξ[1], 𝒫true.ξ[end]])
@rput df
@rput dftrue
@rput BI
R"""
library(tidyverse)
p <- ggplot() +
geom_path(data=df, mapping=aes(x=iterate,y=parameter,colour=type)) +
geom_hline(data=dftrue, mapping = aes(yintercept=parameter,colour=type),linetype="dashed") +
 facet_wrap(~type,scales="free") + theme_light()
png("~/.julia/dev/DataFusion/figs/traceplots.png")
	show(p)
dev.off()
p

print(df %>% filter(iterate>BI) %>% group_by(type) %>%  summarize(m = mean(parameter)))
print(dftrue %>% group_by(type) %>%  summarize(m = mean(parameter)))
"""

θξ = ec(θ,2)
pmξ = [mean(ec(θξ,i)) for i ∈ eachindex(θξ[1])]
println(pmξ)
pmξ1 = θ[end][2]
pmξ2 = θ[end-1][2]
pmξ3 = θ[100][2]

tgr = collect(0:.01:1.0)
pp  = plot(tgr, [μ(tgr[i], pmξ, 𝒫.J) for i ∈ eachindex(tgr)])
plot!(pp, tgr, [μ(tgr[i], pmξ1, 𝒫.J) for i ∈ eachindex(tgr)])
plot!(pp, tgr, [μ(tgr[i], pmξ2, 𝒫.J) for i ∈ eachindex(tgr)])
plot!(pp, tgr, [μ(tgr[i], pmξ3, 𝒫.J) for i ∈ eachindex(tgr)])
png(pp, "/Users/Frank/.julia/dev/DataFusion/figs/averageshape.png")
