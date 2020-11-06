using DataFrames
using Kalman, GaussianDistributions, LinearAlgebra
using GaussianDistributions: ⊕
using Plots
using Distributions
using StaticArrays

cd(@__DIR__)
include("ffbs.jl")

ec(x,i) = map(u->u[i],x)



## Example 1
𝒫 = 0
A(k,𝒫) = [-0.4 0.2; -0.1 0.8]
a(k,𝒫) = [0.0; 0.0]
Q(k,𝒫) = 0.1*[0.2 0.0; 0.0 .2]
H(k,𝒫) = [1.0 0.0]
R(k,𝒫) = 0.01*Matrix(1.0I, 1, 1)
# ------------------ generate data
N = 250 # generate N data, starting from x0 (which is excluded from the dataset)
d = size(A(1,𝒫))[1]
x0 = rand(Gaussian(zeros(2), [0.1 0.0; 0.0 0.1]))
x = [x0]
y =[] # irrelevant but of right type # observations are at indices 2:N
for k in 2:N+1
     global x
     xk = A(k-1,𝒫)*x[k-1] + a(k-1,𝒫) + rand(Gaussian(zeros(d), Q(k-1,𝒫)))
     push!(x,xk)
     yk =  rand(Gaussian(H(k,𝒫)*x[k], R(k,𝒫)))
     push!(y,yk)
end


# Smoothing
d = 2
m0= zeros(d) ; P0=0.0*Matrix(1.0I, d, d)
(m, P), (m⁻, P⁻) = ff(y, (m0,P0),𝒫)
(mˢ, Pˢ) = bsmooth(y, (m, P), (m⁻, P⁻),𝒫)

pp1 = Plots.plot()
Plots.plot!(pp1, 0:N, ec(x,1),label="data",alpha=0.3)
Plots.plot!(pp1, 1:N, ec(m,1),label="filtered")
Plots.plot!(pp1, 1:N, ec(mˢ,1),label="smoothed")
pp2 = Plots.plot()
Plots.plot!(pp2, 0:N, ec(x,2),label="data",alpha=0.3)
Plots.plot!(pp2, 1:N, ec(m,2),label="filtered")
Plots.plot!(pp2, 1:N, ec(mˢ,2),label="smoothed")
plot(pp1,pp2, layout=@layout [a;b])


# ffbs
M = 10 # nr of ffbs pats
Ys = [bsample(y, (m, P), (m⁻, P⁻), 𝒫) for _ ∈ 1:M]
p1 = Plots.plot(); p2 = Plots.plot()
for i in 1:M
    Plots.plot!(p1, 1:N, ec(Ys[i],1),alpha=0.4, label="", color=["steelblue"])
    Plots.plot!(p2, 1:N, ec(Ys[i],2),alpha=0.4, label="", color=["steelblue"])
end
Plots.plot!(p1, 0:N, ec(x,1),markersize=1,label="",colour=:orange)
Plots.plot!(p2, 0:N, ec(x,2),markersize=1,label="",colour=:orange)
plot(p1,p2, layout=@layout [a;b])


## Example 3.6 Sarkka car motion

# Model equations
# x[k] = A(k-1) * x[k-1] + a[k] +  N(0,Q[k-1])
# y[k] = H(k) * x[k] + N(0,R(k))
Δ = 0.01
A(k,𝒫) = [1.0 0.0 Δ 0.0 ;0.0 1.0 0.0 Δ ;0.0 0.0 1.0 0.0 ;0.0 0.0 0.0 1.0 ]
a(k,𝒫) = [0.0; 0.0; 0.0; 0.0]
q1=.2; q2 =0.6
Q(k,𝒫) = [q1*Δ^3/3.0 0.0 q1*Δ^2/2.0 0.0; 0.0 q2*Δ^3/3.0 0.0 q2*Δ^2/2.0; q1*Δ^2/2.0 0.0 q1*Δ 0.0 ; 0.0 q2*Δ^2/2.0 0.0 q2*Δ]
H(k,𝒫) = [1.0 0.0 0.0 0.0 ; 0.0 1.0 0.0 0.0]
R(k,𝒫) = 0.02*Matrix(1.0I, 2, 2)


# ------------------ generate data
N = 1000
d = size(A(1,𝒫))[1]
x0 = zeros(d)
x = [x0]
y =[] # irrelevant but of right type # observations are at indices 2:N
for k in 2:N+1
    println(k)
     global x
     xk = A(k-1,𝒫)*x[k-1] + a(k-1,𝒫) + rand(Gaussian(zeros(d), Q(k-1,𝒫)))
     push!(x,xk)
    yk =  rand(Gaussian(H(k,𝒫)*x[k], R(k,𝒫)))
     push!(y,yk)
end
# Smoothing
m0= zeros(d) ; P0=0.0*Matrix(1.0I, d, d)
(m, P), (m⁻, P⁻) = ff(y, (m0,P0),𝒫)
(mˢ, Pˢ) = bsmooth(y, (m, P), (m⁻, P⁻), 𝒫)

p1 = Plots.plot()
Plots.plot!(p1, 0:N, ec(x,1),label="data",alpha=0.6)
Plots.plot!(p1, 1:N, ec(m,1),label="filtered")
Plots.plot!(p1, 1:N, ec(mˢ,1),label="smoothed")
p2 = Plots.plot()
Plots.plot!(p2, 0:N, ec(x,2),label="data",alpha=0.6)
Plots.plot!(p2, 1:N, ec(m,2),label="filtered")
Plots.plot!(p2, 1:N, ec(mˢ,2),label="smoothed")
p3 = Plots.plot()
Plots.plot!(p3, 0:N, ec(x,3),label="data")
Plots.plot!(p3, 1:N, ec(m,3),label="filtered")
Plots.plot!(p3, 1:N, ec(mˢ,3),label="smoothed")
p4 = Plots.plot()
Plots.plot!(p4, 0:N, ec(x,4),label="data")
Plots.plot!(p4, 1:N, ec(m,4),label="filtered")
Plots.plot!(p4, 1:N, ec(mˢ,4),label="smoothed")

p = Plots.plot(size=(1400, 800))
outp = plot(p1,p2,p3,p4,size=(1400, 800), layout=@layout [a ;b; c; d])
#png(outp, "test_kf.png")



M = 10 # nr of ffbs pats
Ys = [bsample(y, (m, P), (m⁻, P⁻), 𝒫) for _ ∈ 1:M]
p1 = Plots.plot(); p2 = Plots.plot();p3 = Plots.plot();p4 = Plots.plot()
for i in 1:M
    Plots.plot!(p1, 1:N, ec(Ys[i],1),alpha=0.4, label="", color=["steelblue"])
    Plots.plot!(p2, 1:N, ec(Ys[i],2),alpha=0.4, label="", color=["steelblue"])
    Plots.plot!(p3, 1:N, ec(Ys[i],3),alpha=0.4, label="", color=["steelblue"])
    Plots.plot!(p4, 1:N, ec(Ys[i],4),alpha=0.4, label="", color=["steelblue"])
end
Plots.scatter!(p1, 0:N, ec(x,1),markersize=1,label="")
Plots.scatter!(p2, 0:N, ec(x,2),markersize=1,label="")
Plots.scatter!(p3, 0:N, ec(x,3),markersize=1,label="")
Plots.scatter!(p4, 0:N, ec(x,4),markersize=1,label="")
outp = plot(p1,p2,p3,p4, layout=@layout [a b;c d])
png(outp, "test_bffs.png")
