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
This method creates a SPAC node to use in `ForestTower2021` project:

    create_spac(proj::ForestTower2021{FT},
                site::String = "NiwotRidge",
                sm::AbstractStomatalModel = OSMWang{FT}(),
                vcmax::FT = FT(30),
                kmax::FT = FT(0.3),
                chl::FT = FT(40)
    ) where {FT<:AbstractFloat}

Create a SPACMono struct , given
- `proj` [`ForestTower2021`](@ref) type project identifier
- `site` The flux tower site identifier, must be `NiwotRidge` or `Ozark`
- `sm` `AbstractStomatalModel` type stomatal model. Default is `OSMWang`.
- `vcmax` Maximum carboxylation rate as an input. Default is `30`.
- `kmax` Maximum whole plant hydraulic conductance as an input. Default is
    `0.3`.
- `chl` Chlorophyll content. Default is `40`.
"""
create_spac(proj::ForestTower2021{FT},
            site::String = "NiwotRidge",
            sm::AbstractStomatalModel = OSMWang{FT}(),
            vcmax::FT = FT(30),
            kmax::FT = FT(0.3),
            chl::FT = FT(40)
) where {FT<:AbstractFloat} =
(
    # make sure that site is supported
    @assert site in ["NiwotRidge", "Ozark"];

    if site == "NiwotRidge"
        _soil_hs = VanGenuchten{FT}(stype = "Cambisol",
                                        α = 129.66,
                                        n = 1.35301,
                                       Θs = 0.477,
                                       Θr = 0.196);
        _tree_hs = create_tree(FT(-1), FT(6), FT(12.5), FT[0,-0.5,-1],
                               collect(FT,0:0.5:13));
        for _root in _tree_hs.roots
            _root.sh = deepcopy(_soil_hs);
        end;

        _node = SPACMono{FT}(soil_bounds = FT[0,-0.5,-1],
                              air_bounds = collect(FT,0:0.5:13),
                                  z_root = -1,
                                z_canopy = 12.5,
                                plant_hs = _tree_hs,
                                      ga = 32.088,
                                      la = 128.352,
                                latitude = 40.03,
                               longitude = -105.55,
                               elevation = 3050,
                           stomata_model = sm);
        for _iPS in _node.plant_ps
            _iPS.g_min   = 0.001;
            _iPS.g_min25 = 0.001;
            _iPS.g_max   = 0.08;
            _iPS.g_max25 = 0.08;
        end;
        _node.canopy_rt.Ω       = 0.48;
        _node.canopy_rt.clump_a = 0.48;
        update_LAI!(_node, FT(4));
        update_Weibull!(_node, FT(4.09), FT(5.82));
    elseif site=="Ozark"
        # soil_hs = create_soil_VC(VanGenuchten{FT}(), "Silt Loam");
        _soil_hs = VanGenuchten{FT}(stype = "Ozark",
                                        α = 1.368,
                                        n = 2.6257,
                                       Θs = 0.45,
                                       Θr = 0.067);
        _tree_hs = create_tree(FT(-1), FT(9), FT(18.5), FT[0,-0.5,-1],
                               collect(FT,0:0.5:19));
        for _root in _tree_hs.roots
            _root.sh = deepcopy(_soil_hs);
        end;

        _node = SPACMono{FT}(soil_bounds = FT[0,-0.5,-1],
                              air_bounds = collect(FT,0:0.5:19),
                                  z_root = -1,
                                z_canopy = 18.5,
                                plant_hs = _tree_hs,
                                      ga = 413.223,
                                      la = 1735.537,
                                latitude = 38.74,
                               longitude = -92.2,
                               elevation = 219.4,
                           stomata_model = sm);
        for _iPS in _node.plant_ps
            _iPS.g_min   = 0.001;
            _iPS.g_min25 = 0.001;
            _iPS.g_max   = 0.08;
            _iPS.g_max25 = 0.08;
        end;
        _node.canopy_rt.Ω       = 0.69;
        _node.canopy_rt.clump_a = 0.69;
        update_LAI!(_node, FT(4.2));
        update_Weibull!(_node, FT(5.703), FT(0.953));
    end;

    # update fitted Vcmax, Kmax, and Chlrophyll content
    update_VJRWW!(_node, vcmax);
    update_Kmax!(_node, kmax);
    update_Cab!(_node, chl);
    initialize_spac_canopy!(_node);

    return _node
)
