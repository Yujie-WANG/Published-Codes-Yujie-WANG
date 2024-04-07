module ResearchProjects

using LazyArtifacts

using CanopyLayers: SIF_740, SIF_fluxes!, canopy_fluxes!, canopy_geometry!, canopy_matrices!, fluspect!, short_wave!
using DataFrames: DataFrame
using Photosynthesis: AbstractPhotoModelParaSet, AirLayer, C3CLM, GCO₂Mode, Leaf, PCO₂Mode, Q10TD, leaf_ETR!, leaf_fluorescence!, leaf_photosynthesis!, photo_TD_from_set
using PkgUtility: K_STEFAN, M_H₂O, lower_quadratic, nanmean, nanstd, numerical∫, parse_timestamp, read_csv, tinfo
using PlantHydraulics: LeafHydraulics, SteadyStateMode, VanGenuchten, pressure_profile!, soil_p_25_swc, temperature_effects!
using PlotPlants: create_canvas, latex_symbol, latex_unit, save_canvas!, set_titles!, set_xlabels!, set_xticks!, set_xticklabels!, set_xylabels!, set_xylims!, set_xyticks!, set_ylabels!, set_ylims!,
      use_serif_tex
using ProgressMeter: @showprogress
using PyPlot: figure
using SoilPlantAirContinuum: SPACMono, initialize_spac_canopy!, update_Cab!, update_Kmax!, update_LAI!, update_VJRWW!, update_Weibull!, zenith_angle
using Statistics: mean
using StomataModels: AbstractStomatalModel, CanopyLayer, EmpiricalStomatalModel, GswDrive, OSMWang, gas_exchange!, gsw_control!, stomatal_conductance, update_leaf_TP!
using UnPack: @unpack
using WaterPhysics: latent_heat_vapor, saturation_vapor_pressure


# export 2021_canopy_complexity types and functions
export CanopyComplexity2021, plot_CC2021!


# General files
include("general/constants.jl"  )
include("general/projects.jl"   )
include("general/simulations.jl")
include("general/create_spac.jl")
include("general/query_data.jl" )

# Project CanopyComplexity2021
include("2021_canopy_complexity/diurnal.jl")
include("2021_canopy_complexity/figures.jl")
include("2021_canopy_complexity/spectra.jl")


end # module
