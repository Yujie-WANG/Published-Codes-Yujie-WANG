###############################################################################
#
# Create SPAC Simple or Mono
#
###############################################################################
"""
    create_spac(proj::NocturnalGS2020{FT},
                ele::FT = FT(1400)
    ) where {FT<:AbstractFloat}

Create SPACSimple or SPACMono struct , given
- `proj` [`LeafInvestment2020`](@ref) or [`FluxTower2020`](@ref) type project
- `ele` Elevation. Optional, default at Flagstaff
"""
function create_spac(
            proj::NocturnalGS2020{FT},
            ele::FT = FT(1400)
) where {FT<:AbstractFloat}
    # create node
    _node = SPACSimple{FT}(
                        gaba      = 1000  ,
                        laba      = 4758.5,
                        lai       = 4.7585,
                        elevation = ele   ,
                        width     = 0.1   );
    _node.hs.leaf = LeafHydraulics{FT}(
                        vc=WeibullSingle{FT}(b=1.897, c=2.203),
                        k_sla=0.0176);
    _node.hs.root = RootHydraulics{FT}(
                        vc=WeibullSingle{FT}(b=1.879, c=2.396),
                        k_max=34.383);
    _node.hs.stem = StemHydraulics{FT}(
                        vc=WeibullSingle{FT}(b=2.238, c=9.380),
                        k_max=76.029);
    _node.ps.Vcmax25 = 61.74 ;
    _node.ps.Jmax25  = 111.13;
    _node.ps.Rd25    = 1.5;

    return _node
end
