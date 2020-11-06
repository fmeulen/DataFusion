# Set observation scheme
recording = (
	P = Periodic1D(θinit_full...),  # still don't understand why this has to be passed to the function (also θtrue is unknown normally)
	obs = load_data(
		ObsScheme(
			LinearGsnObs(
				0.0, (@SVector [0.0]);
				L=(@SMatrix [1.0]), Σ=(@SMatrix [0.01])
			)
		),
		data
	),
	t0 = 0.0,
	#x0_prior = GsnStartingPt([0.0], [4.0^2]),#
	x0_prior = KnownStartingPt(y1),
)

if false
# Now with two data sources

# type 1
η = 0.01
t, v = 1.0, [1.0]
obs1 = LinearGsnObs(t, v; L=[1.0], Σ = reshape(η*[1.0], (1,1)),full_obs=true)

# type 2, i.e. another observation scheme:
t, v = 1.0, [1.0; 1.0]
obs2 = LinearGsnObs(t, v; L=[1.0 0.0; 0.0 1.0], Σ = η*[1.0 0.0; 0.0 1.0],full_obs=true)

# data defined externally:
tt, xx = [1.0, 2.0, 4.0], [[1.], [2.5 ;3.5], [10.0] ]
# template of an observation scheme in pattern: obs, obs2, obs, obs2, obs, ...
obs_scheme = ObsScheme(obs1, obs2; pattern=[1,2,1])
# decorate the data
observs = load_data(obs_scheme, tt, xx)

recording = (
	P = Periodic1D(θinit_full...),  # still don't understand why this has to be passed to the function (also θtrue is unknown normally)
	obs = observs,
	t0 = 0.0,
	x0_prior = GsnStartingPt([0.0], [4.0^2]),#KnownStartingPt(y1),
)

end
