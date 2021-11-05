module ResearchProjects

using LazyArtifacts
using Revise

using BenchmarkTools: @benchmark
using CanopyLayers: EVI, EVI2, LSWI, LeafBios, NDVI, NIRv, SIF_740, SIF_757,
      SIF_771, SIF_fluxes!, WaveLengths, canopy_fluxes!, canopy_geometry!,
      canopy_matrices!, create_canopy_rt, create_leaf_bios, create_wave_length,
      create_rt_dims, dladgen, e2phot, fluspect!, initialize_rt_module,
      short_wave!
using ConstrainedRootSolvers: BisectionMethod, ReduceStepMethod,
      ReduceStepMethodND, SolutionTolerance, SolutionToleranceND, find_peak,
      find_zero
using DataFrames: DataFrame
using Distributed: @everywhere, addprocs, pmap, rmprocs, workers
using GriddingMachine: AbstractDataset, CanopyHeightGLAS, ClumpingIndexMODIS,
      ERA5LandHourly, GPPMPIv006, GPPVPMv20, ERA5SingleLevelsHourly,
      GriddedDataset, LAIMODISv006, LandElevation, LandMaskERA5,
      LeafChlorophyll, LeafSLAButler, SIFTropomi740, TreeDensity,
      VGMAlphaJules, VGMLogNJules, VGMThetaRJules, VGMThetaSJules, WoodDensity,
      fetch_ERA5!, lat_ind, load_LUT, lon_ind, mask_LUT!, query_LUT, read_LUT,
      regrid_LUT, save_LUT!
using NCDatasets: Dataset
using Photosynthesis: AbstractPhotoModelParaSet, AirLayer, C3CLM, C4CLM,
      GCO₂Mode, Q10TD, leaf_photosynthesis!, photo_TD_from_set
using PkgUtility: GRAVITY, K_STEFAN, M_H₂O, nanmax, nanmean, nanmin,
      nanpercentile, nanstd, numerical∫, parse_timestamp, read_csv, read_nc,
      rmse, save_csv!, save_nc!, send_email!, terror, tinfo, twarn, ρ_H₂O,
      ρg_MPa
using PlantHydraulics: BrooksCorey, GrassLikeOrganism, LeafHydraulics,
      LogisticSingle, NonSteadyStateMode, PalmLikeOrganism, PowerSingle,
      SteadyStateMode, StemHydraulics, TreeLikeOrganism, VanGenuchten,
      WeibullSingle, create_grass, create_palm, create_soil_VC, create_tree,
      critical_flow, end_pressure, pressure_profile!, roots_flow!,
      soil_k_ratio_swc, soil_p_25_swc, temperature_effects!, update_PVF!,
      vc_integral, xylem_k_ratio, xylem_risk
using PlotPlants: calculate_density, create_canvas, latex_symbol, latex_unit,
      line_regress, plot_density!, plot_hexbin!, plot_line_regress!,
      preview_data, save_canvas!, save_gif!, set_titles!, set_xlabels!,
      set_xlims!, set_xticklabels!, set_xticks!, set_xylabels!, set_xylims!,
      set_xyticklabels!, set_xyticks!, set_ylabels!, set_ylims!,
      set_yticklabels!, test_slope, use_serif_tex
using ProgressMeter: @showprogress
using PyPlot: figure
using SoilPlantAirContinuum: SPACMono, initialize_spac_canopy!, layer_fluxes!,
      update_Cab!, update_LAI!, update_Kmax!, update_VJR!, update_VJRWW!,
      update_Weibull!, zenith_angle
using StomataModels: AbstractStomatalModel, CanopyLayer, ESMBallBerry,
      ESMMedlyn, EmpiricalStomatalModel, GswDrive, OSMWang, gas_exchange!,
      gsw_control!, stomatal_conductance, update_leaf_TP!
using UnPack: @unpack
using WaterPhysics: latent_heat_vapor, relative_surface_tension,
      relative_viscosity, saturation_vapor_pressure




# export general functions
export dynamic_workers!, reprocess_FLUXNET!

# export 2021_forest_tower types and functions
export ForestTower2021, divide_results!, fit_tower!, plot_FT2021!,
       plot_FT2021_SI!, sif_series!




# General files
include("general/constants.jl"   )
include("general/projects.jl"    )
include("general/flux_tower.jl"  )
include("general/threadings.jl"  )
include("general/create_spac.jl" )

# Project ForestTower2021
include("2021_forest_tower/figures.jl"   )
include("2021_forest_tower/fit_traits.jl")
include("2021_forest_tower/match_data.jl")
include("2021_forest_tower/query_data.jl")
include("2021_forest_tower/simulation.jl")
include("2021_forest_tower/test_patm.jl" )




end # module
