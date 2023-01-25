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

export ReactionRates, PlasmaDistributions, Integrations, Testing

include("reaction-rates.jl")

"""
    Read eqdsk.
    Define domain R,Z,
    Build ψ(r,z)
    Read excel and build n(ψ) and T(ψ).
    Compute I(r,z).
"""
module PlasmaDistributions
    using DataFrames, Interpolations, XLSX
    using EQDSKReader
    using ..ReactionRates

    export Content, DDN, DT, σv
    export AbstractDistribution, DDDistribution, DTDistribution
    export domain, ψ, n, σv, Ti, I
    export concentrations
    export load_excel, load_psi

    """
        dR_dV(n, σv)

    Reactivity density, DD case, cm-3s-1
    """
    dR_dV(n, σv) = @. 0.5 * n^2 * σv

    """
        dR_dV(n, σv)

    Reactivity density, DT case, cm-3s-1
    """
    dR_dV(n1, n2, σv) = @. n1 * n2 * σv

    """
    Base type for DD(pure D) and DT(DT mix) distributions.
    """
    abstract type AbstractDistribution end

    domain(d::AbstractDistribution) = d.rmin, d.rmax, d.zmin, d.zmax
    ψ(d::AbstractDistribution) = d.ψ
    n(d::AbstractDistribution) = d.n
    n(d::AbstractDistribution, ψ) = n(d)(ψ)
    n(d::AbstractDistribution, r, z) = n(d, ψ(d)(r, z))
    Ti(d::AbstractDistribution) = d.T
    Ti(d::AbstractDistribution, ψ) = Ti(d)(ψ)
    Ti(d::AbstractDistribution, r, z) = Ti(d, ψ(d)(r, z))
    σv(d::AbstractDistribution) = d.σv
    σv(d::AbstractDistribution, ψ) = σv(d)(Ti(d)(ψ))
    σv(d::AbstractDistribution, r, z) = σv(d, ψ(d)(r, z))

    # Be careful with broadcasting over r, z given as matrices.
    # Either the two matrices should be of the same shape (like meshgrid)
    # or represent a grid points.

    """
        I(d::AbstractDistribution) -> I(ψ)

    ## Returns
        - function of neutron source intensity vs. ψ
    """
    I(d::AbstractDistribution) = ψ -> dR_dV(n(d, ψ), σv(d, ψ))  # DD like rate by default

    """
        I(d::AbstractDistribution, ψ) -> Float64 or Array{Float64}

    ## Returns
        - neutron source intensity, cm^-3⋅c^-1,
          either scalar or array of the same shape as ψ
    """
    I(d::AbstractDistribution, ψ) = I(d)(ψ)

    I(d::AbstractDistribution, r, z) = I(d, ψ(d)(r, z))

    struct DDDistribution <: AbstractDistribution
        rmin::Float64
        rmax::Float64
        zmin::Float64
        zmax::Float64
        ψ::Function
        T::Function
        n::Function
        σv::Function
    end

    function DDDistribution(c::Content, df::DataFrame)
        DDDistribution(
            load_psi(c)..., create_t_vs_psi(df), create_n_vs_psi(df), ReactionRates.σv(DDN)
        )
    end

    struct DTDistribution <: AbstractDistribution
        rmin::Float64
        rmax::Float64
        zmin::Float64
        zmax::Float64
        ψ::Function
        T::Function
        n::Function
        σv::Function
        fuel_ratio::Float64
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

    # """
    # Function of ψ fuel ratio.
    # """
    # function _do_get_concentrations(ψ, n, fuel_ratio)
    #     _n = n(ψ)
    #     nt = fuel_ratio.(ψ) .* _n
    #     nd = _n - nt
    #     nd, nt
    # end

    concentrations(d::DTDistribution) = ψ -> _do_get_concentrations(ψ, d.n, d.fuel_ratio)
    concentrations(d::DTDistribution, ψ) = concentrations(d)(ψ)
    concentrations(d::DTDistribution, r, z) = concentrations(d)(ψ(d)(r, z))

    """
        I(d::DTDistribution) -> Function

    Specialization to compute source intensity for DT case, cm-3s-1

    ## Returns
        Function to compute neutron source strengs as functions of ψ, cm^-3s-1.
    """
    function I(d::DTDistribution)
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

using .PlasmaDistributions

module Integrations

    using Cuba
    using ..PlasmaDistributions

    export total_yield, torroidal_segment_yield, integrate_torroidal_segment
    export torroidal_segment_moment_0, torroidal_segment_moment_1

    """
    total_yield(d::AbstractDistribution) -> Float64

    Compute total neutron yield for distribution `d`.

    Calculates integral of I(r,z) in cylinder coordinates (R,Z,ϕ),
    where
    - R - major radius
    - Z - height
    - ϕ - torroidal angle

    ## Arguments
    - d - plasma distribution specification

    ## Returns
    - total neutron yield, s^-1
    - error, s^-1
    - number of estimation
    - fail or not fail
    """
    function total_yield(d::AbstractDistribution)
        rmin, rmax, zmin, zmax = domain(d)
        integrand(r, z) = I(d, r, z)
        integrate_torroidal_segment(integrand, (rmin, rmax), (zmin, zmax))
    end

    # function torroidal_segment_yield(I, rbin, zbin)
    #     integral, error, neval, fail = integrate_torroidal_segment(I, rbin, zbin)
    #     # Normalization:
    #     # - 1e6   - m^3/cm^3
    #     integral * 1e6, error * 1e6, neval, fail
    # end

    function integrate_torroidal_segment(I, rbin, zbin)
        R0 = rbin[1]
        ΔR = rbin[2] - R0
        Z0 = zbin[1]
        ΔZ = zbin[2] - Z0
        function integrand(x, f)
            r = ΔR * x[1] + R0
            z = ΔZ * x[2] + Z0
            f[1] = r * I(r, z)
        end
        integral, error, _, neval, fail, _ = cuhre(integrand)
        # Normalization:
        # - 2π    - integral over torroidal angle,
        # - 1e6   - m^3/cm^3
        # - ΔR⋅ΔZ - area of R,Z integration domain, square meters
        normalization = 2e6π * ΔR * ΔZ
        integral[1] * normalization, error[1] * normalization, neval, fail
    end

    torroidal_segment_moment_0(I, rbin, zbin) = integrate_torroidal_segment(I, rbin, zbin)

    function torroidal_segment_moment_1(I, rbin, zbin)
        R0 = rbin[1]
        ΔR = rbin[2] - R0
        Z0 = zbin[1]
        ΔZ = zbin[2] - Z0
        function integrand(x, f)
            r = ΔR * x[1] + R0
            z = ΔZ * x[2] + Z0
            f[:] = [r, z] .* (r * I(r, z))
        end
        integral0, error0, neval, fail = integrate_torroidal_segment(I, rbin, zbin)
        fail == 0 || raise(ArgumentError("Cannot integrate with default args"))
        integral, error, _, neval, fail, _ = cuhre(integrand, 2, 2)
        # Normalization:
        # - 2π    - integral over torroidal angle,
        # - 1e6   - m^3/cm^3
        # - ΔR⋅ΔZ - area of R,Z integration domain, square meters
        normalization = 2e6π * ΔR * ΔZ / integral0
        integral .* normalization, error .* normalization .+ error0 / integral0, neval, fail
    end

end

using .Integrations

module Testing

    using ..PlasmaDistributions

    export TestDistribution

    """
        struct TestDistribution

    Distribution for testing and code optimization.
    Represents cylinder ring with 1.5 ≤ r ≤ 2.5, -0.5 ≤ z ≤ 0.5.
    ψ(r,z) is 1.0 if (r,z) ∈ the ring, 0 - otherwise.
    Total yield is the volume of the ring in cm^3: 1e6*2π*(2.5^2 - 1.5^2)/2.
    """
    struct TestDistribution <: AbstractDistribution end

    PlasmaDistributions.domain(::TestDistribution) = 1.5, 2.5, -0.5, 0.5
    PlasmaDistributions.ψ(::TestDistribution) = (r, z) -> in_domain.(r, z)
    PlasmaDistributions.n(::TestDistribution) = ψ -> 1.0 .- ψ
    PlasmaDistributions.Ti(::TestDistribution) = ψ -> 1.0 .- ψ
    PlasmaDistributions.I(::TestDistribution) = ψ -> 1.0 .- ψ
    in_domain(r, z) = (1.5 <= r <= 2.5 && -0.5 <= z <= 0.5 ? 0.0 : 1.0)

end

end
