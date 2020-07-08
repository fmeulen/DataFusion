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
p2=Plots.plot( ec(θ,3),label="c")
p3=Plots.plot( ec(θ,4),label="σ2")
p4=Plots.plot( first.(ec(θ,5)),label="ψ1")
p5=Plots.plot( last.(ec(θ,5)),label="ψ2")
Plots.plot(ec(θ,1)./ec(θ,4), label="α/σ2")
p = Plots,plot(p1,p2,p3,p4,p5)
show(p)
#png(p, "pars.png")


df = DataFrame(iterate= repeat(1:ITER,5),
	parameter= vcat(ec(θ,1),ec(θ,3),ec(θ,4),first.(ec(θ,5)),last.(ec(θ,5))),
	type=repeat(["alpha","c","sigma2","psi1","psi2"],inner=ITER))
dftrue = DataFrame(type=["alpha","c","sigma2","psi1","psi2"], parameter=[𝒫true.α, 𝒫true.c, 𝒫true.σ2, 𝒫true.ψ[1], 𝒫true.ψ[2]])
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
