

BI = div(ITER,2)

postmean_paths = ec1([mean(map(x->x[i],X[BI:ITER-1])) for i in eachindex(X[1])])
writedlm("../../csv/postmean_paths.csv", postmean_paths)

θξ = ec(θ,2)[BI:ITER]
pmξ = [mean(ec(θξ,i)) for i ∈ eachindex(θξ[1])]
tgr = collect(0:.001:1.0)
dfmupost = DataFrame(t=tgr, y=[μ(tgr[i], pmξ, 𝒫.J) for i ∈ eachindex(tgr)])


df = DataFrame(iterate= repeat(1:ITER,6),
	parameter= vcat(ec(θ,1),ec(θ,3),first.(ec(θ,4)),last.(ec(θ,4)),first.(ec(θ,2)),last.(ec(θ,2))),
	type=repeat(["alpha","sigma2","psi1","psi2","xi1","xilast"],inner=ITER))
dftrue = DataFrame(type=["alpha","sigma2","psi1","psi2","xi1","xilast"], parameter=[𝒫true.α, 𝒫true.σ2, 𝒫true.ψ[1], 𝒫true.ψ[2], 𝒫true.ξ[1], 𝒫true.ξ[end]])
dfpath = DataFrame(t=t[2:end], postmean=postmean_paths, y=ec1(y))
@rput df
@rput dftrue
@rput BI
@rput dfpath
@rput dfmupost
R"""
library(tidyverse)
p <- ggplot() +
geom_path(data=df, mapping=aes(x=iterate,y=parameter,colour=type)) +
geom_hline(data=dftrue, mapping = aes(yintercept=parameter,colour=type),linetype="dashed") +
 facet_wrap(~type,scales="free") + theme_light() + theme(legend.position='none')
pdf("~/.julia/dev/DataFusion/figs/traceplots.pdf",width=7,height=4)
	show(p)
dev.off()
p

ppath <- dfpath %>% ggplot() +
geom_point(aes(x=t,y=y),colour='blue',size=0.3,alpha=0.6)+
geom_path(aes(x=t,y=postmean),size=0.3) +
 ylab("concentration") + xlab("time elapased") + theme_light()
pdf("~/.julia/dev/DataFusion/figs/posterior_path.pdf",width=7,height=4)
	show(ppath)
dev.off()

pshape <- dfmupost %>% ggplot(aes(x=t,y=y)) + geom_path()+ theme_light()
pdf("~/.julia/dev/DataFusion/figs/mupost.pdf",width=7,height=4)
	show(pshape)
dev.off()


print(df %>% filter(iterate>BI) %>% group_by(type) %>%  summarize(m = mean(parameter)))
print(dftrue %>% group_by(type) %>%  summarize(m = mean(parameter)))
"""

# informal jl plotting

if false
	pl = Plots.scatter(𝒫.t[2:end],first.(y), markersize=3)
	for k ∈ eachindex(X[end-10:end])
	    Plots.plot!(pl,𝒫.t[2:end],ec1(X[k]),label="")
	end
	Plots.plot!(pl, 𝒫.t[2:end],postmean_paths)
	display(pl)
	#png(pl,"test1.png")
end




p1=Plots.plot( ec(θ,1),label="α")
p2=Plots.plot( first.(ec(θ,2)),label="ξ1")
p3=Plots.plot( ec(θ,3),label="σ2")
p4=Plots.plot( first.(ec(θ,4)),label="ψ1")
p5=Plots.plot( last.(ec(θ,4)),label="ψ2")
p6=Plots.plot( ec(θ,1)./ec(θ,3),label="α/σ2")
println(𝒫true.α/𝒫true.σ2)
