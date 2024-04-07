###############################################################################
#
# Modify SPAC Simple or Mono
#
###############################################################################
"""
Modify `SPACMono` to meet specific requirements.
"""
function mod_spac! end




"""
This method sets a `SPACMono` to default settings without initialize it:
- LAI = 3
- Kmax = 10
- ZSA = VZA = RAZ = 0
- Clear sky with transmittance of 0.7

    mod_spac!(node::SPACMono{FT}) where {FT<:AbstractFloat}

Modify the SPAC structure, given
- `node` A `SPACMono` type structure
"""
mod_spac!(node::SPACMono{FT}, sza::FT=FT(0)) where {FT<:AbstractFloat} =
(
    # update the angles
    node.angles.sza = sza;
    node.angles.vza = 0;
    node.angles.raa = 0;

    # update direct and diffuse radiation based on SZA
    _in_rad_bak = node.in_rad;
    _e_all_dire = sum(_in_rad_bak.E_direct  .* node.wl_set.dWL) / 1000;
    _e_all_diff = sum(_in_rad_bak.E_diffuse .* node.wl_set.dWL) / 1000;
    _sb         = 1367 * 0.7 ^ (1/cosd(sza)) * cosd(sza);
    _sd         = 1367 * 0.3 * (1 - 0.7 ^ (1/cosd(sza))) * cosd(sza);
    node.in_rad.E_direct  .*= _sb / _e_all_dire;
    node.in_rad.E_diffuse .*= _sd / _e_all_diff;

    return nothing
)




"""
This method applies a vertical Vcmax and Cab gradient for a SPACMono

    mod_spac!(node::SPACMono{FT}, c3c4::String)

Changes the Vcmax and Cab for each canopy layers, given
- `node` A `SPACMono` type structure
- `c3c4` A indicator of C3 plant or C4 plant
"""
mod_spac!(node::SPACMono{FT}, c3c4::String) where {FT<:AbstractFloat} = (
    _lai = node.la / node.ga;
    for _i_can in eachindex(node.plant_ps)
        _x_rt = exp( -0.15 * _lai * (1 - _i_can / length(node.plant_ps)) );
        node.plant_ps[_i_can].ps.Vcmax25 *= _x_rt;
        node.plant_ps[_i_can].ps.Jmax25  *= _x_rt;
        node.plant_ps[_i_can].ps.Rd25    *= _x_rt;
        node.plant_ps[_i_can].ps.Vcmax   *= _x_rt;
        node.plant_ps[_i_can].ps.Jmax    *= _x_rt;
        node.plant_ps[_i_can].ps.Rd      *= _x_rt;
        # ! these impacts APAR, do not change them this way
        #node.leaves_rt[_i_can].Cab      *= _x_rt;
        #node.leaves_rt[_i_can].Car      *= _x_rt;
        fluspect!(node.leaves_rt[_i_can], node.wl_set);
    end;

    return nothing
)








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
This method create a default SPACMono and initialize it.

    create_spac(FT)

Create a SPACMono structure, given
- `FT` Floating number type
"""
create_spac(FT, sm::AbstractStomatalModel=OSMWang{FT}()) =
(
    _node = SPACMono{FT}(soil_bounds = FT[0,-0.5,-1],
                          air_bounds = collect(FT,0:0.1:3),
                              z_root = -0.4,
                            z_canopy = 1.99,
                       stomata_model = sm);
    mod_spac!(_node);

    # initialize the node
    initialize_spac_canopy!(_node);

    return _node
)




"""
This method creates a SPAC node with different number of leaves and layers for
    [`CanopyComplexity2021`](@ref) project:

    create_spac(proj::CanopyComplexity2021{FT}) where {FT<:AbstractFloat}

Create a `SPACMono` structure, given
- `proj` [`CanopyComplexity2021`](@ref) type structure
"""
create_spac(proj::CanopyComplexity2021{FT},
            mode::String,
            sm::AbstractStomatalModel = OSMWang{FT}()
) where {FT<:AbstractFloat} =
(
    # create SPAC and initialize it
    if mode == "IJKX"
        _node = create_spac(FT, sm);
    elseif mode == "2KX"
        _node = create_spac(FT, sm);
        _node.plant_ps = [CanopyLayer{FT}(n_leaf=2) for _i in 1:20];
        mod_spac!(_node);
    elseif mode == "KX"
        _node = create_spac(FT, sm);
        _node.plant_ps = [CanopyLayer{FT}(n_leaf=1) for _i in 1:20];
        mod_spac!(_node);
    elseif mode == "2X"
        _pss  = [CanopyLayer{FT}(n_leaf=1) for _i in 1:2];
        _node = SPACMono{FT}(soil_bounds = FT[0,-0.5,-1],
                              air_bounds = collect(FT,0:0.1:3),
                                  z_root = -0.4,
                                z_canopy = 0.19,
                                plant_ps = _pss,
                           stomata_model = sm);
        mod_spac!(_node);
    elseif mode == "1X"
        _pss  = [CanopyLayer{FT}(n_leaf=1)];
        _node = SPACMono{FT}(soil_bounds = FT[0,-0.5,-1],
                              air_bounds = collect(FT,0:0.1:3),
                                  z_root = -0.4,
                                z_canopy = 0.09,
                                plant_ps = _pss,
                           stomata_model = sm);
        mod_spac!(_node);
    end;

    return _node
)




"""
This method creates a SPAC node with different number of leaves and layers for
    [`CanopyComplexity2021`](@ref) project:

    create_spac(proj::CanopyComplexity2021{FT}) where {FT<:AbstractFloat}

Create a `SPACMono` structure, given
- `proj` [`CanopyComplexity2021`](@ref) type structure
"""
create_spac(proj::CanopyComplexity2021{FT},
            sm::AbstractStomatalModel = OSMWang{FT}(),
            gradients::Bool = true,
            hour::Int = 12
) where {FT<:AbstractFloat} =
(
    # create nodes to work on
    _node_ijkx = create_spac(proj, "IJKX");
    _node_2kx  = create_spac(proj, "2KX");
    _node_kx   = create_spac(proj, "KX");
    _node_2x   = create_spac(proj, "2X");
    _node_1x   = create_spac(proj, "1X");
    _nodes     = [_node_ijkx, _node_2kx, _node_kx, _node_2x, _node_1x];
    _soil_hs   = VanGenuchten{FT}(stype = "Ozark", α = 1.368, n = 2.6257,
                                  Θs = 0.45, Θr = 0.067);
    for _node in _nodes
        for _root in _node.plant_hs.roots
            _root.sh = deepcopy(_soil_hs);
        end;
        _node.ga = 413.223;
        _node.la = 1735.537;
        _node.latitude = 38.74;
        _node.longitude = -92.2;
        _node.elevation = 219.4;
        _node.canopy_rt.Ω = 0.69;
        _node.canopy_rt.clump_a = 0.69;
        for _iPS in _node.plant_ps
            _iPS.g_min   = 0.001;
            _iPS.g_min25 = 0.001;
            _iPS.g_max   = 0.095;
            _iPS.g_max25 = 0.095;
        end;
        update_LAI!(_node, FT(4.2));
        update_VJRWW!(_node, FT(75));
        update_Weibull!(_node, FT(5.703), FT(0.953));
        update_Kmax!(_node, FT(1.55));
        update_Cab!(_node, FT(57.23));
        if gradients
            mod_spac!(_node_ijkx, "C3");
            mod_spac!(_node_2kx, "C3");
            mod_spac!(_node_kx, "C3");
        end;
    end;

    # use Ozark data for the sensitivity analysis
    _rawdf = query_data(proj, 2019, "Ozark");
    _mask  = (177 .<= _rawdf.Day .< 178) .*
             (_rawdf.Hour .== hour) .*
             (0 .<= _rawdf.Minu .< 30);
    _newdf = _rawdf[_mask,:];

    # update soil water matrices for all nodes
    _w_soil = FT(_newdf.SWC[1]) / 100;
    _ψ_soil = soil_p_25_swc(_node_ijkx.plant_hs.roots[1].sh, _w_soil);
    for _node in _nodes
        for _root in _node.plant_hs.roots
            _root.p_ups = _ψ_soil;
        end;
    end;

    # update PAR related information
    _zenith = zenith_angle(_node_ijkx.latitude, FT(_newdf.Day[1]),
                           FT(_newdf.Hour[1]), FT(_newdf.Minu[1]));
    _zenith = min(88, _zenith);
    for _node in _nodes
        mod_spac!(_node, _zenith);
        (_node).angles.sza = _zenith;
    end;

    # use Ozark radiation
    _in_rad_bak = deepcopy(_node_ijkx.in_rad);
    _in_SW = (numerical∫(_in_rad_bak.E_direct, _node_ijkx.wl_set.dWL) +
              numerical∫(_in_rad_bak.E_diffuse, _node_ijkx.wl_set.dWL)) /
              1000;
    _ratio = _newdf.IN_RAD[1] ./ _in_SW;
    for _node in _nodes
        _node.in_rad.E_direct  .= _in_rad_bak.E_direct  .* _ratio;
        _node.in_rad.E_diffuse .= _in_rad_bak.E_diffuse .* _ratio;
    end;

    # create an environment struct to use in all modes
    _envir   = AirLayer{FT}();
    _envir = deepcopy(_node_ijkx.envirs[1]);
    _envir.t_air = _newdf.T_AIR[1] + 273.15;
    _envir.p_atm = _newdf.P_ATM[1] * 1000;
    _envir.p_a   = _envir.p_atm * 4e-4;
    _envir.p_O₂  = _envir.p_atm * 0.209;
    _envir.p_sat = saturation_vapor_pressure(_envir.t_air);
    _envir.vpd   = _newdf.VPD[1] * 100;
    _envir.p_H₂O = _envir.p_sat - _envir.vpd;
    _envir.RH    = _envir.p_H₂O / _envir.p_sat;
    _envir.wind  = _newdf.WIND[1];

    # prescribe leaf temperature
    _tl = (_newdf.LW_OUT[1] / 0.97 / K_STEFAN() ) ^ 0.25;
    for _node in _nodes
        for _iPS in _node.plant_ps
            _iPS.T = _tl;
        end;
    end;

    return _nodes, _envir
)
