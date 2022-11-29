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

include("reaction-rates.jl")

module Reactivity

    dR_dV(n, σv) = 0.5 * n^2 * σv

    dR_dV(n1, n2, σv) = n1 * n2 * σv

end

"""
    Read eqdsk.
    Define domain R,Z,
    Build ψ(r,z)
    Read excel and build n(ψ) and T(ψ).
    Compute I(r,z).
"""
module PlasmaDistribution
    using EQDSKReader, DataFrames, Interpolations, XLSX
    using ..ReactionRates: DDN, DT, σv
    using ..Reactivity: dR_dV

    export Content, DDN, DT, σv
    export Distribution, load_psi, load_excel, create_t_vs_psi, T, n, Idd, Idt

    function load_psi(c::Content)
        rmin, rmax = extrema(c.rbbbs)
        zmin, zmax = extrema(c.zbbbs)
        rmin, rmax, zmin, zmax, create_normalized_psi_interpolator(c)
    end

    load_psi(path::AbstractString) = load_psi(Content(path))

    function load_excel(path)::DataFrame
        data = Float64.(XLSX.readdata(path, 1, "A5:D183")[:, [1, 3, 4]])
        d = DataFrame(data, ["ψ", "T", "n"])
        pushfirst!(d, [0.0, d[1, Not(1)]...])  # Add row for ψ == 0.0
        select!(d, :ψ, :T, :n => (x -> 1e13x) => :n) # scale concentation to cm-3
    end

    function create_interpolator(df::DataFrame, s::Symbol)
        interpolate(df.ψ, df[!, s], SteffenMonotonicInterpolation())
    end

    create_t_vs_psi(df::DataFrame) = create_interpolator(df, :T)
    create_n_vs_psi(df::DataFrame) = create_interpolator(df, :n)

    struct Distribution
        rmin
        rmax
        zmin
        zmax
        ψ
        T
        n
        σvdd
        σvdt
    end

    function Distribution(c::Content, df::DataFrame)
        Distribution(
            load_psi(c)..., create_t_vs_psi(df), create_n_vs_psi(df), σv(DDN), σv(DT)
        )
    end
    """
        T(d::Distribution, r, z)

    2D temperature distibution, keV

    # Arguments
    - d - structure specifying plasma distribution
    - r - the location major radius(radiuses)
    - z - ... height (heights)
    """
    T(d::Distribution, r, z) = d.T(d.ψ(r, z))

    """
        n(d::Distribution, r, z)

    2D ion concentration distibution, cm-3

    # Arguments
    - d - structure specifying plasma distribution
    - r - the location major radius(radiuses)
    - z - ... height (heights)
    """
    n(d::Distribution, r, z) = d.n(d.ψ(r, z))

    """
        default_compute_fuel_ratio(r) -> Function

    Create default fuel ratio computing method - just return initial constant

    # Returns
        - Function of ψ to compute fuel ratio in dependency of ψ
    """
    default_compute_fuel_ratio(ratio) = ψ -> ratio

    """
        Idd(d::Distribution, ψ)

    DD plasma neutron source intensity over ψ, cm-3s-1.

    # Arguments
    - d - structure specifying plasma distribution
    - r - the location major radius(radiuses)
    - z - ... height (heights)
    """
    function Idd(d::Distribution, ψ)
        1.0 < ψ && return 0.0
        dR_dV(d.n(ψ), d.σvdd(d.T(ψ)))
    end

    """
        Idd(d::Distribution, r, z)

    2D DD plasma neutron source intensity, cm-3s-1.

    # Arguments
    - d - structure specifying plasma distribution
    - r - the location major radius(radiuses)
    - z - ... height (heights)
    """
    function Idd(d::Distribution, r, z)
        !(d.rmin <= r <= d.rmax && d.zmin <= z <= d.zmax) && return 0.0
        Idd(d, d.ψ(r, z))
    end

    """
        Idt(d::Distribution, ψ, η::Function=default_compute_fuel_ratio(0.5))

    DT plasma neutron source intensity over ψ, cm-3s-1.

    # Arguments
    - d - structure specifying plasma distribution
    - r - the location major radius(radiuses)
    - z - ... height (heights)
    - η - function specifying fuel ratio distribution (default, const 0.5)
    """
    function Idt(d::Distribution, ψ, η::Function=default_compute_fuel_ratio(0.5))
        1.0 < ψ && return 0.0
        n = d.n(ψ)
        ratio = η(ψ)
        n_d, n_t = (1.0 - ratio)n, ratio * n
        dR_dV(n_d, n_t, d.σvdt(d.T(ψ)))
    end

    """
        Idt(d::Distribution, r, z, η::Function=default_compute_fuel_ratio(0.5))

    2D DT plasma neutron source intensity, cm-3s-1.

    # Arguments
    - d - structure specifying plasma distribution
    - r - the location major radius(radiuses)
    - z - ... height (heights)
    - η - function specifying fuel ratio distribution
    """
    function Idt(d::Distribution, r, z, η::Function=default_compute_fuel_ratio(0.5))
        !(d.rmin <= r <= d.rmax && d.zmin <= z <= d.zmax) && return 0.0
        Idt(d, d.ψ(r, z), η)
    end

end

end
