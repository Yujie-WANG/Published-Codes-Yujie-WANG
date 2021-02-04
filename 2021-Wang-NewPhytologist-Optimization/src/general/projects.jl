###############################################################################
#
# Abstract type of projects
#
###############################################################################
"""
    abstract type AbstractProject{FT}

Hierarchy of AbstractProject:
- [`CanopyWater2020`](@ref)
- [`ClumpingFactor2020`](@ref)
- [`FluxTower2020`](@ref)
- [`GriddingMachine2021`](@ref)
- [`LeafInvestment2020`](@ref)
- [`LeafInvestment2021`](@ref)
- [`NocturnalGS2020`](@ref)
- [`OptimalGM2020`](@ref)
- [`SIFComparison2020`](@ref)
"""
abstract type AbstractProject{FT} end




"""
    struct NocturnalGS2020{FT}
"""
struct NocturnalGS2020{FT} <: AbstractProject{FT} end
