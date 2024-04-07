module ResearchProjects

using CanopyLayers: LeafBios, WaveLengths, canopy_fluxes!, canopy_geometry!,
      canopy_matrices!, create_canopy_rt, create_leaf_bios, create_rt_dims,
      create_wave_length, fluspect!, short_wave!, SIF_740, SIF_757, SIF_771,
      SIF_fluxes!
using ConstrainedRootSolvers: ReduceStepMethodND, SolutionToleranceND,
      find_peak
using DataFrames: DataFrame
using PkgUtility: nanmean, numerical∫, read_nc, save_csv!, tinfo
using PlotPlants: create_canvas, set_titles!, set_xlabels!, set_ylabels!,
      use_serif_tex
using SoilPlantAirContinuum: SPACMono, initialize_spac_canopy!, update_LAI!,
      zenith_angle
using UnPack: @unpack




# export 2020_clumping_factor types and functions
export ClumpingFactor2020, fit_leafbio, sif_simulation!




# General files
include("general/projects.jl"   )
include("general/create_spac.jl")

# Project ClumpingFactor2020
include("2020_clumping_factor/fit_leafbio.jl" )
include("2020_clumping_factor/simulate_sif.jl")




end # module
