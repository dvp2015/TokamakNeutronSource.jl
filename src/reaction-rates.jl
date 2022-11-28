"""
Module implements reaction rates parametrization algorithm of Bosch&Hale.

## References
> [1] _Improved formulas for fusion cross-sections and thermal particles_
>     H.-S. Bosch and G.M. Hale, Nuclear Fusion vol. *32*, No 4 (1992), p624
>     https://iopscience.iop.org/article/10.1088/0029-5515/32/4/I07
"""
module ReactionRates

export DDN, DT, σv, σv_ddn, σv_dt, sigmav_ddn, sigmav_dt

"""
    Traits to select reaction rate algorithms.
"""
abstract type CollidedParticles end
"""
    Trait to select T<d,n>4He reaction rate algorithm.
"""
struct DT <: CollidedParticles end
"""
    Trait to select D<d,n>3He reaction rate algorithm.
"""
struct DDN <: CollidedParticles end

"""
    σv(::Type{<:CollidedParticles}) -> Function

Generic implementation.

# Arguments
- tp - trait to select DT reaction rate algorithm
"""
function σv(tp::Type{<:CollidedParticles})
    error("σv is not implemented for type $tp")
end

"""
    σv_dt(t) -> Float::64

Implementation of T(d,n)4He reaction rate, cm3/s

# Arguments
- t - temperature, keV


"""
function σv_dt(t::Real)::Float64
    b_g_sq = 34.3827^2
    m_r_c = 1124656
    c1 = 1.17302e-9
    c2 = 1.51361e-2
    c3 = 7.51886e-2
    c4 = 4.60643e-3
    c5 = 1.35000e-2
    c6 = -1.06750e-4
    c7 = 1.36600e-5
    t_min = 0.2

    t < t_min && return 0.0
    θ = t / (1 - (t * (c2 + t * (c4 + t * c6))) / (1 + t * (c3 + t * (c5 + t * c7))))
    ξ = ∛(b_g_sq / (4.0θ))
    c1 * θ * sqrt(ξ / (m_r_c * t^3)) * exp(-3.0ξ)
end

"""
    σv(::Type{DT}) -> Function

T(d,n)4He reaction rate, cm3/s

# Arguments
- DT - trait to select DT reaction rate algorithm


"""
σv(::Type{DT}) = σv_dt

"""
    sigmav_dt(t) -> Float64

T(d,n)4He reaction rate, cm3/s, alias


# Arguments
- t - temperature, keV


"""
sigmav_dt = σv_dt

"""
    σv_ddn(t) -> Float::64

Implementation of D(d,n)3He reaction rate, cm3/s

# Arguments
- t - temperature, keV


"""
function σv_ddn(t)::Float64
    b_g_sq = 31.3970^2
    m_r_c = 937814
    c1 = 5.43360e-12
    c2 = 5.85778e-3
    c3 = 7.68222e-3
    c5 = -2.96400e-6
    t_min = 0.2

    t < t_min && return 0.0
    θ = t / (1 - t * c2 / (1 + t * (c3 + t * c5)))
    ξ = ∛(b_g_sq / (4.0θ))
    c1 * θ * sqrt(ξ / (m_r_c * t^3)) * exp(-3.0ξ)
end

"""
    σv(::Type{DDN}) -> Function

D(d,n)3He reaction rate, cm3/s

# Arguments
- DDN - trait to select DDN reaction rate algorithm


"""
σv(::Type{DDN}) = σv_ddn

"""
sigmav_ddn(t) -> Float::64

D(d,n)3He reaction rate, cm3/s, alias

# Arguments
- t - temperature, keV


"""
sigmav_ddn = σv_ddn

# TODO dvp: implement reaction rates for 3He(d,p)4He and D(d,p)T reactions.
# This is necessary to estimate all the 'ordinary thermal' neutron sources.
end
