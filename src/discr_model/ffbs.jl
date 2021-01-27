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
function ff(y, (m0,P0), 𝒫)
    Tm = typeof(m0)
    TP = typeof(P0)
    m, P = Tm[], TP[]
    m⁻, P⁻ = Tm[], TP[]

    mprev = m0; Pprev = P0
    for k ∈ eachindex(y)
        #println(k)
        push!(m⁻, A(k-1,𝒫) * mprev + a(k-1,𝒫))
        push!(P⁻, A(k-1,𝒫) * Pprev * A(k-1,𝒫)' + Q(k-1,𝒫))
        v = y[k] - (H(k,𝒫) * m⁻[k])
        S = H(k,𝒫) * P⁻[k] * H(k,𝒫)' + R(k,𝒫)
        K = (P⁻[k] * H(k,𝒫)')/S
        mprev = m⁻[k] + K * v
        Pprev = P⁻[k] - K * S * K'
        push!(m, mprev)
        push!(P, Pprev)
    end
    (m, P), (m⁻, P⁻)
end


function ff!(y, (m0,P0), (m, P), (m⁻, P⁻), 𝒫)
    mprev = m0; Pprev = P0
    for k ∈ eachindex(y)
        m⁻[k] =  A(k-1,𝒫) * mprev + a(k-1,𝒫)
        P⁻[k] = A(k-1,𝒫) * Pprev * A(k-1,𝒫)' + Q(k-1,𝒫)
        v = y[k] .- (H(k,𝒫) * m⁻[k])
        S = H(k,𝒫) * P⁻[k] * H(k,𝒫)' + R(k,𝒫)
        K = (P⁻[k] * H(k,𝒫)')/S
        m[k] = mprev = m⁻[k] + K * v
        P[k] = Pprev = P⁻[k] - K * S * K'
    end
    nothing
end


"""
    bsample(y, (m, P), (m⁻, P⁻))

Backward sampling to generate a sample path from the smoothing distribution
y: data
(m, P): filtered mean and covariances
(m⁻,P⁻): forward predicted mean and covariances (compute in forward filtering)

Returns: a sample path from the smoothing distribution
"""
function bsample((m, P), (m⁻, P⁻), 𝒫)
    #dump(typeof(m))
    yout = 0.0*m  #zeros(typeof(m),length(m))
    bsample!(yout,(m, P), (m⁻, P⁻), 𝒫)
    yout
end

"""
    bsample!(yout,(m, P), (m⁻, P⁻), 𝒫)

Inplace version of bsample
"""
function bsample!(yout,(m, P), (m⁻, P⁻), 𝒫)
    # FIXME: the additive drift term a(k) should maybe also come in here
    N = length(m)
    yout[N] = rand(Gaussian(m[N],P[N]))
    for k ∈ N-1:-1:1
        G = (P[k] * A(k,𝒫)')/P⁻[k+1]
        z = m[k] + G * (yout[k+1] - m⁻[k+1])
        cv = P[k] - G *  P⁻[k+1] * G'
        yout[k] = rand(Gaussian(z,cv))
    end
    nothing
end
