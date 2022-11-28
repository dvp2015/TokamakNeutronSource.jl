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
    using EQDSKReader, XLSX

    function load_psi(path)
        c = Content(path)
        rmin, rmax = extrema(c.rbbbs)
        zmin, zmax = extrema(c.zbbbs)
        (rmin=rmin, rmax=rmax, zmin=zmin, zmax=zmax, ψ=create_normalized_psi_interpoloator(c))
    end

# function

end

end
