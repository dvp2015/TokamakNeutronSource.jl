"""
    TokamakNeutronSource

Code to compute neutron emission intensity for given plasma state.

## Description
The plasma state is presented with:
- `ψ(r,z)` - magnetic surface coordinate spatial distribution
- `n(ψ)` - ion concentration in dependency of magnetic surface coordinate `ψ`
- `T(ψ)` - Tritium concentration in dependency of magnetic surface coordinate `ψ`
- `r` - ratio of atomic concentraions of tritium to deuterium one

"""
module TokamakNeutronSource

"""
Module implements reaction rates parametrization algorithm of Bosch&Hale.

## References
> [1] _Improved formulas for fusion cross-sections and thermal particles_
>     H.-S. Bosch and G.M. Hale, Nuclear Fusion vol. *32*, No 4 (1992), p624
>     https://iopscience.iop.org/article/10.1088/0029-5515/32/4/I07
"""
module ReactionRates

    export σv_ddn, σv_dt, sigmav_ddn, sigmav_dt

    abstract type AbstractParameters end

    function _θ(t::Float64, p::AbstractParameters)::Float64 end

    @inline function σv(t, p::AbstractParameters)::Float64
        t < p.t_min && return 0.0
        __θ = _θ(t, p)
        ξ = ∛(p.b_g_sq / (4.0__θ))
        p.c1 * __θ * sqrt(ξ / (p.m_r_c * t^3)) * exp(-3ξ)
    end

    struct ParametersDT <: AbstractParameters
        b_g_sq::Float64
        m_r_c::Float64
        c1::Float64
        c2::Float64
        c3::Float64
        c4::Float64
        c5::Float64
        c6::Float64
        c7::Float64
        t_min::Float64

        function ParametersDT()
            new(
                34.3827^2,
                1124656,
                1.17302e-9,
                1.51361e-2,
                7.51886e-2,
                4.60643e-3,
                1.35000e-2,
                -1.06750e-4,
                1.36600e-5,
                0.2,
            )
        end
    end

    const DT = ParametersDT()

    _θ(t::Float64, p::ParametersDT)::Float64 =
        t / (
            1 -
            (t * (p.c2 + t * (p.c4 + t * p.c6))) / (1 + t * (p.c3 + t * (p.c5 + t * p.c7)))
        )

    """
        σv_dt(t) -> Float::64

    \$T(d,n)^{4}He\$ reaction rate, \$cm^{3}/s\$

    # Arguments
    - t - temperature, keV


    """
    σv_dt(t) = σv(t, DT)
    sigmav_dt = σv_dt
    struct ParametersDDN <: AbstractParameters
        b_g_sq::Float64
        m_r_c::Float64
        c1::Float64
        c2::Float64
        c3::Float64
        c5::Float64
        t_min::Float64

        function ParametersDDN()
            new(
                31.3970^2,
                937814,
                5.43360e-12,
                5.85778e-3,
                7.68222e-3,
                -2.96400e-6,
                0.2,
            )
        end
    end

    function σv_dt(tarray::AbstractArray{T,N}) where {T<:Real,N}
        σv_dt.(tarray)
    end

    const DDN = ParametersDDN()

    _θ(t::Float64, p::ParametersDDN)::Float64 =
        t / (1 - t * p.c2 / (1 + t * (p.c3 + t * p.c5)))

    """
        σv_ddn(t) -> Float::64

    \$D(d,n)^{3}He\$ reaction rate, \$cm^{3}/s\$

    # Arguments
    - t - temperature, keV


    """
    σv_ddn(t) = σv(t, DDN)
    sigmav_ddn = σv_ddn

    function σv_ddn(tarray::AbstractArray{T,N}) where {T<:Real,N}
        σv_ddn.(tarray)
    end

# TODO dvp: implement reaction rates for 3He(d,p)4He and D(d,p)T reactions
end

end
