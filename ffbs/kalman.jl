using DataFrames
using RCall
using Kalman, GaussianDistributions, LinearAlgebra
using GaussianDistributions: ⊕
using Plots
meancov(g) = g.μ, g.Σ
symmetrize(g) = Gaussian(g.μ, Symmetric(g.Σ) + 0.0001*I)

# forward simulate
x0 = [1., 0.]
P0 = Matrix(1.0I, 2, 2)

#Φ = [0.8 0.5; -0.1 0.8]
Φ = [0.8 0.2; -0.1 0.8]
b = zeros(2)
Q = [0.2 0.0; 0.0 1.0]

ymean = [0.0]
H = [1.0 0.0]
R = Matrix(1.0I, 1, 1)

# sample

N = 10
x = rand(Gaussian(x0, P0))
xs = [x]
ys = Any[]
for i in 1:N
     global x
     x = Φ*x + b + rand(Gaussian(zero(x), Q))
     y = H*x + rand(Gaussian(ymean, R))
     push!(xs, x)
     push!(ys, y)
end


if false


#--- check with M's implementation
Φ = A(1)
b = a(1)

# filter
P0=5.0*Matrix(1.0I, 2, 2)
p = Gaussian(x0, P0)
ps = [p]
ppreds = [p]
for i in 1:N
    global p, ps
    p = Φ*p + b
    p = p ⊕ Gaussian(zero(x0), Q(1))
    push!(ppreds, p)

    p, yres, S = Kalman.correct(Kalman.JosephForm(), p, (Gaussian(y[i], R(1)), H(1)))
    push!(ps, p)
    #correct(method::JosephForm, u::T, (v, H)::Tuple{<:Gaussian, <:Any}) where T
end


function backward_rand_kernel(x, Gf, Gpred, Phi, b) where T
    xf, Pf = meancov(Gf)
    xpred, Ppred = meancov(Gpred)

    J = Pf*Phi'/Ppred
    xs = xf +  J*(x - (Phi*xf  + b))
    Ps = Pf - J*Ppred*J' # as Ppred = M.Phi*P*M.Phi' + M.Q
    println(Ps)
    rand(symmetrize(Gaussian(xs, Ps)))
end


x = rand(symmetrize(ps[end]))
X = [x]
for i in N-1:-1:1
    println("i $i, $x")
    global x
    x = backward_rand_kernel(x, ps[i], ppreds[i], Φ, b)
    push!(X, x)
end


end


## own implementation

# x[k] = A(k-1) * x[k-1] + a[k] +  N(0,Q[k-1])
# y[k] = H(k) * x[k] + N(0,R(k))

A(k) = [-0.4 0.2; -0.1 0.8]
a(k) = [0.0; 0.0]
Q(k) = 0.1*[0.2 0.0; 0.0 .2]
H(k) = [1.0 0.0]
R(k) = 0.01*Matrix(1.0I, 1, 1)

 # H(k) = [1.0 0.0; 1.0 1.0]
 # R(k) = 0.1*Matrix(1.0I, 2, 2)

d = size(A(1))[1]

# prior mean and variance
m = [[0.0; 0.0]]
P = [5.0*Matrix(1.0I, 2, 2)]

# generate data
N = 500
x0 = rand(Gaussian(m[1], P[1]))
x = [x0]
y =[H(1)*x0] # irrelevant but of right type # observations are at indices 2:N

for k in 2:N
     global x
     xk = A(k-1)*x[k-1] + a(k-1) + rand(Gaussian(zeros(d), Q(k-1)))
     push!(x,xk)
     yk =  rand(Gaussian(H(k)*x[k], R(k)))
     push!(y,yk)
end


function ffbs(y, (m0,P0)) # first el of y is irrelevant, first index is for the prior
    m, P = [m0], [P0]
    m⁻, P⁻ = [m0], [P0]
    N = length(y)
    # KF
    for k ∈ 2:N
        push!(m⁻, A(k-1) * m[k-1] + a(k-1))
        push!(P⁻, A(k-1) * P[k-1] * A(k-1)' + Q(k-1))

        v = y[k] - H(k) * m⁻[k]
        S = H(k) * P⁻[k] * H(k)' + R(k)
        K = P⁻[k] * H(k)'/S
        push!(m, m⁻[k] + K * v)
        push!(P, P⁻[k] - K * S * K')
    end
    # Backw sampling
    mˢ = [m[N]]
    Pˢ = [P[N]]
    for k ∈ N-1:-1:1
        G = (P[k] * A(k)')/P⁻[k+1]
        pushfirst!(mˢ, m[k] + G * (first(mˢ) - m⁻[k+1]))
        pushfirst!(Pˢ, P[k] + G * (first(Pˢ) - P⁻[k+1]) * G')
    end
    (m, P), (mˢ, Pˢ)
end

m0= [0.0; 0.0] ; P0=2.0*Matrix(1.0I, 2, 2)
(m, P), (mˢ, Pˢ) = ffbs(y,(m0,P0))
hcat(x,m, mˢ)

ec(x,i) = map(u->u[i], x)
df  = DataFrame(t= 1:N, x1=ec(x,1), x2 = ec(x,2), f1=ec(m,1), f2=ec(m,2), s1 = ec(mˢ,1), s2 = ec(mˢ,2))
dfout = df[2:end,:]

@rput dfout
R"""
library(tidyverse); library(gridExtra)
g1 <- dfout %>% ggplot() + geom_path(aes(x=t,y=x1))+
    geom_path(aes(x=t,y=f1),colour='red')+
    geom_path(aes(x=t,y=s1),colour='blue')

g2 <-   dfout %>% ggplot() +  geom_path(aes(x=t,y=x2))+
        geom_path(aes(x=t,y=f2),colour='red')+
        geom_path(aes(x=t,y=s2),colour='blue')

grid.arrange(g1,g2)
"""
