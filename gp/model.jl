@diffusion_process Periodic1D{T} begin
    :dimensions
    process --> 1
    wiener --> 1
    :parameters
    (ᾱ, ω̄, ϕ̄, δ, σ̄) --> T
    :additional
    constdiff --> true
end

@diffusion_process Periodic1DAux{K, R} begin
    :dimensions
    process --> 1
    wiener --> 1
    :parameters
    (ᾱ, ω̄, ϕ̄, δ, σ̄) --> K
    :auxiliary_info
    t0 --> Float64
    T --> Float64
    vT --> R

    :additional
    constdiff --> true
    linear --> true
end

ν(t, P::Union{Periodic1D, Periodic1DAux}) = P.δ*pdf(Beta(exp(P.ω̄),exp(P.ϕ̄)),mod(t,1.0))

DD.b(t, x, P::Periodic1D)= @SVector [exp(P.ᾱ)*(ν(t,P) - x[1])]
DD.σ(t, x, P::Periodic1D) = @SMatrix [exp(P.σ̄)]

DD.B(t, P::Periodic1DAux) =  @SMatrix [-exp(P.ᾱ)]
DD.β(t, P::Periodic1DAux) = @SVector [exp(P.ᾱ) * ν(t,P)]
DD.σ(t, P::Periodic1DAux) = @SMatrix [exp(P.σ̄)]

# Declare which parameters are not changing
DD.const_parameter_names(::Type{<:Periodic1D}) = ()
DD.const_parameter_names(::Type{<:Periodic1DAux}) = (:t0, :T, :vT, :xT)
