###############################################################################
#
# Create SPAC Simple or Mono
#
###############################################################################
"""
Create soil-plant-air continuum node.
"""
function create_spac end




"""
This method creates a SPAC node to use in `ClumpingFactor2020` project:

    ccreate_spac(
                proj::ClumpingFactor2020{FT},
                ci::FT
    ) where {FT<:AbstractFloat}

Create a SPACMono struct , given
- `proj` [`ClumpingFactor2020`](@ref) type project identifier
- `ci` Clumping index
"""
create_spac(proj::ClumpingFactor2020{FT}, ci::FT) where {FT<:AbstractFloat} =
(
    _node = SPACMono{FT}(ga = 32.088,
                         la = 128.352,
                   latitude = 40.03,
                  longitude = -105.55,
                  elevation = 3050);
    _node.canopy_rt.Ω       = ci;
    _node.canopy_rt.clump_a = ci;
    update_LAI!(_node, FT(4));
    for _leaf in _node.leaves_rt
        _leaf.N = 1.6;
        _leaf.Cab = 30.0;
        _leaf.Car = 5.0;
        _leaf.Ant = 2.75;
        _leaf.Cs = 0.0;
        _leaf.Cw = 5e-3;
        _leaf.Cm = 0.0;
        _leaf.Cx = 0.0;
        fluspect!(_leaf, _node.wl_set);
    end;
    initialize_spac_canopy!(_node);

    return _node
)
