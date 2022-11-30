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

export ReactionRates, PlasmaDistribution

include("reaction-rates.jl")

"""
    Read eqdsk.
    Define domain R,Z,
    Build ψ(r,z)
    Read excel and build n(ψ) and T(ψ).
    Compute I(r,z).
"""
module PlasmaDistribution
    using EQDSKReader, DataFrames, Interpolations, XLSX
    using ..ReactionRates

    export Content, DDN, DT, σv
    export AbstractDistribution, DDDistribution, DTDistribution
    export Ti, n, σv, I
    export concentrations
    export load_excel, load_psi

    function dR_dV(n, σv)
        @. 0.5 * n^2 * σv
    end

    function dR_dV(n1, n2, σv)
        @. n1 * n2 * σv
    end

    abstract type AbstractDistribution end

    ψ(d::AbstractDistribution) = d.ψ
    n(d::AbstractDistribution) = d.n
    n(d::AbstractDistribution, ψ) = n(d)(ψ)
    n(d::AbstractDistribution, r, z) = n(d, ψ(d)(r, z))
    concentrations(d::AbstractDistribution) = ψ ->  (n(d,ψ),)
    concentrations(d::AbstractDistribution, ψ) = concentrations(d)(ψ)
    concentrations(d::AbstractDistribution, r, z) = concentrations(d, ψ(d)(r,z))
    Ti(d::AbstractDistribution) = d.T
    Ti(d::AbstractDistribution, ψ) = Ti(d)(ψ)
    Ti(d::AbstractDistribution, r, z) = Ti(d, ψ(d)(r, z))
    σv(d::AbstractDistribution) = d.σv
    σv(d::AbstractDistribution, ψ) = σv(d)(Ti(d)(ψ))
    σv(d::AbstractDistribution, r, z) = σv(d, ψ(d)(r, z))
    I(d::AbstractDistribution) = ψ -> dR_dV(n(d, ψ), σv(d, ψ))  # mono particle rate
    I(d::AbstractDistribution, ψ) = I(d)(ψ)
    I(d::AbstractDistribution, r, z) = I(d, ψ(d)(r, z))

    struct DDDistribution <: AbstractDistribution
        rmin
        rmax
        zmin
        zmax
        ψ
        T
        n
        σv
    end

    function DDDistribution(c::Content, df::DataFrame)
        DDDistribution(
            load_psi(c)..., create_t_vs_psi(df), create_n_vs_psi(df), ReactionRates.σv(DDN)
        )
    end

    struct DTDistribution <: AbstractDistribution
        rmin
        rmax
        zmin
        zmax
        ψ
        T
        n
        σv
        fuel_ratio
    end

    function DTDistribution(c::Content, df::DataFrame, fuel_ratio=0.5)
        DTDistribution(
            load_psi(c)...,
            create_t_vs_psi(df),
            create_n_vs_psi(df),
            ReactionRates.σv(DT),
            fuel_ratio,
        )
    end

    """
    Constant fuel ratio.
    """
    function _do_get_concentrations(ψ, n, fuel_ratio::T) where {T<:Real}
        _n = n(ψ)
        nt = fuel_ratio .* _n
        nd = _n - nt
        nd, nt
    end

    """
    Function of ψ fuel ratio.
    """
    function _do_get_concentrations(ψ, n, fuel_ratio)
        _n = n(ψ)
        nt = fuel_ratio.(ψ) .* _n
        nd = _n - nt
        nd, nt
    end

    concentrations(d::DTDistribution) = ψ -> _do_get_concentrations(ψ, d.n, d.fuel_ratio)
    concentrations(d::DTDistribution, ψ) = concentrations(d)(ψ)
    concentrations(d::DTDistribution, r, z) = concentrations(d)(ψ(d)(r, z))

    """
        I(d::DTDistribution) -> Function

    Specialization to compute source intensity for DT case, cm-3s-1

    ## Returns
        Function to compute neutron source strengs as functions of ψ, cm^-3s-1.
    """
    function I(d::T) where {T<:DTDistribution}
        cf = concentrations(d)
        ψ -> dR_dV(cf(ψ)..., σv(d, ψ))
    end

    function load_psi(c::Content)
        rmin, rmax = extrema(c.rbbbs)
        zmin, zmax = extrema(c.zbbbs)
        rmin, rmax, zmin, zmax, create_normalized_psi_interpolator(c)
    end

    load_psi(path::AbstractString) = load_psi(Content(path))

    function load_excel(path)::DataFrame
        data = Float64.(XLSX.readdata(path, 1, "A5:D183")[:, [1, 3, 4]])
        d = DataFrame(data, ["ψ", "T", "n"])
        pushfirst!(d, [0.0, d[1, Not(:ψ)]...])  # Add row for ψ == 0.0
        select!(d, :ψ, :T, :n => (x -> 1e13x) => :n) # scale concentation to cm-3
    end

    function create_interpolator(df::DataFrame, s::Symbol)
        itp = extrapolate(interpolate(df.ψ, df[!, s], SteffenMonotonicInterpolation()), 0.0)
        ψ -> itp.(ψ)
    end

    create_t_vs_psi(df::DataFrame) = create_interpolator(df, :T)
    create_n_vs_psi(df::DataFrame) = create_interpolator(df, :n)

end

using .PlasmaDistribution

end
