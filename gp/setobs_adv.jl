v = @SVector [0.0]
t = 0.0
L = @SMatrix [1.0]
obs1 = LinearGsnObs(t, v; L = L, θ=[σ1], Tag=1)

v = @SVector [0.0]
t = 0.0
L = @SMatrix [1.0]
obs2 = LinearGsnObs(t, v; L = L, θ=[σ2], Tag=2)

v = @SVector [0.0, 0.0]
t = 0.0
L = @SMatrix [1.0; 1.0]
obs3 = LinearGsnObs(t, v; L = L, θ=[σ1, σ2], Tag=3)

OBS.Σ(o::LinearGsnObs{1}) = SDiagonal{1,Float64}((o.θ[1]^2,)) # NOTE the dispatch on the `tag`s value
OBS.Σ(o::LinearGsnObs{2}) = SDiagonal{1,Float64}((o.θ[1]^2,))
OBS.Σ(o::LinearGsnObs{3}) = SDiagonal{2,Float64}((o.θ[1]^2, o.θ[2]^2))

# data defined externally:
dat_all = CSV.read("observations.csv")


dat = dat_all[1:K,:]

#tt, xx = [1.0, 2.0, 4.0], [[1.], [2.5 ;3.5], [10.0] ]
# obs_scheme = ObsScheme(obs1, obs2. obs3; pattern=[1,2,1])


X = Vector{Float64}[]
for r in eachrow(dat)
	if r[:obsscheme]=="obs1"
		u = tryparse(Float64,r[:chl_water])
		push!(X,[u])
	elseif r[:obsscheme]=="obs2"
		u = tryparse(Float64,r[:chl])
		push!(X,[u])
	elseif r[:obsscheme]=="obs3"
		u = tryparse.(Float64,[r[:chl_water], r[:chl]])
		push!(X,u)
	end
end
obss = String.(dat[!,:obsscheme])
pat = (obss.=="obs1") + 2*(obss.=="obs2") + 3*(obss.=="obs3")
obs_scheme = ObsScheme(obs1,obs2,obs3; pattern=pat)
ts = Float64.(dat[!,:time_elapsed])


observs = load_data(obs_scheme, ts, X)
#summary.(observs)

y1 = @SVector [X[1][1]]
recording = (
	P = Periodic1D(θinit_full...),  # still don't understand why this has to be passed to the function (also θtrue is unknown normally)
	obs = observs,
	t0 = 0.0,
	x0_prior = KnownStartingPt(y1)#GsnStartingPt([0.0], [4.0^2]),#KnownStartingPt(y1),
)




# Now, if you want to update parameters, then you should simply change all relevant θ values in each relevant observations.
# The relevant observations are part of `GuidProp` structs, so you basically need to iterate through your `GuidProp` structs `PP`
# say and update the `obs` field accordingly, i.e. change values `PP[i].obs.θ[j]` for relevant `GuidProp` struct `i` and parameter `j`.
