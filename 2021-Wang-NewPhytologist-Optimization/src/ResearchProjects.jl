module ResearchProjects




using CLIMAParameters
using CLIMAParameters.Planet
using ConstrainedRootSolvers
using CSV
using DataFrames
using LazyArtifacts
using Parameters
using Photosynthesis
using Pkg.Artifacts
using PlantHydraulics
using PlotPlants
using PyPlot
using SoilPlantAirContinuum
using StomataModels
using WaterPhysics




# constants
struct EarthParameterSet <: AbstractEarthParameterSet end
const EARTH               = EarthParameterSet();
K_STEFAN(FT=Float64)      = FT( Stefan() );
MOLMASS_WATER(FT=Float64) = FT( molmass_water(EARTH) );




# export project types
export AbstractProject,
       NocturnalGS2020

# export general functions
export create_spac

# export 2020_nocturnal_gs functions
export optimal_en_gn,
       plot_gswn_vs_time,
       plot_model_comparison,
       plot_model_extension,
       plot_model_framework,
       plot_model_prediction,
       plot_si_t_leaf




# General files
include("general/constants.jl"  )
include("general/projects.jl"   )
include("general/create_spac.jl")




# project NocturnalGS2020
include("2020_nocturnal_gs/figures.jl"  )
include("2020_nocturnal_gs/gain_risk.jl")




end # module
