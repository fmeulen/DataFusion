# Model equations
# x[k] = A(k-1) * x[k-1] + a[k-1] +  N(0,Q[k-1])
# y[k] = H(k) * x[k] + N(0,R(k))

# See pages 57 and 136 in Sarkka: Bayesian filtering and smoothing
# https://users.aalto.fi/~ssarkka/pub/sde_book.pdf
# (m0,P0): prior mean and cov for the state

"""
    ff(y, (m0,P0))

Forward Kalman filtering
y: data
(m0,P0): prior mean and covariance of state

Returns: (m, P), (m⁻, P⁻), where (m,P) are (mean,covmatrix) filtering distribution
"""
function ff(y, (m0,P0), 𝒫) # forward filter
    m, P = [], []
    m⁻, P⁻ = [], []
    N = length(y)
    # Forw filtering
    for k ∈ 1:N
        if k==1
            push!(m⁻, A(k-1,𝒫) * m0 + a(k-1,𝒫))
            push!(P⁻, A(k-1,𝒫) * P0 * A(k-1,𝒫)' + Q(k-1,𝒫))
        else
            push!(m⁻, A(k-1,𝒫) * m[k-1] + a(k-1,𝒫))
            push!(P⁻, A(k-1,𝒫) * P[k-1] * A(k-1,𝒫)' + Q(k-1,𝒫))
        end
        v = y[k] - H(k,𝒫) * m⁻[k]
        S = H(k,𝒫) * P⁻[k] * H(k,𝒫)' + R(k,𝒫)
        K = (P⁻[k] * H(k,𝒫)')/S
        push!(m, m⁻[k] + K * v)
        push!(P, P⁻[k] - K * S * K')
    end
    (m, P), (m⁻, P⁻)
end


"""
    bsmooth(y, (m0,P0), (m, P), (m⁻, P⁻)) # backw smoothing

Backward Kalman smoothing
y: data
(m, P): filtered mean and covariances
(m⁻,P⁻): forward predicted mean and covariances (compute in forward filtering)

Returns: (mˢ, Pˢ) which are (mean,covmatrix) smoothing distributions
"""
function bsmooth(y, (m, P), (m⁻, P⁻),𝒫) # backw smoothing
    # FIXME: the additive drift term a(k) should maybe also come in here
    N = length(m)
    mˢ = [m[N]]
    Pˢ = [P[N]]
    for k ∈ N-1:-1:1
        G = (P[k] * A(k,𝒫)')/P⁻[k+1]
        pushfirst!(mˢ, m[k] + G * (first(mˢ) - m⁻[k+1]))
        pushfirst!(Pˢ, P[k] + G * (first(Pˢ) - P⁻[k+1]) * G')
    end
    (mˢ, Pˢ)
end

"""
    bsample(y, (m, P), (m⁻, P⁻)) # backw sampling

Backward sampling to generate a sample path from the smoothing distribution
y: data
(m, P): filtered mean and covariances
(m⁻,P⁻): forward predicted mean and covariances (compute in forward filtering)

Returns: a sample path from the smoothing distribution
"""
function bsample(y, (m, P), (m⁻, P⁻), 𝒫) # backw sampling
    # FIXME: the additive drift term a(k) should maybe also come in here
    N = length(m)
    yout = [rand(Gaussian(m[N], Symmetric(P[N])))]
    for k ∈ N-1:-1:1
        G = (P[k] * A(k,𝒫)')/P⁻[k+1]
        z = m[k] + G * (first(yout) - m⁻[k+1])
        cv = Symmetric(P[k] - G *  P⁻[k+1] * G' )
        pushfirst!(yout, rand(Gaussian(z,cv)))
    end
    yout
end
