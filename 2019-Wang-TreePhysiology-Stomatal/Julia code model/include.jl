# this file is used to include every julia script to use to the function

#include the files in

# the math modules
include("./math/weibull_k.jl")
include("./math/weibull_plc.jl")
include("./math/weibull_p_crit.jl")

# the physics modules
include("./physics/energy_boundary_layer_conductance.jl")
include("./physics/energy_emittance_blackbody.jl")
include("./physics/energy_emittance_constants.jl")
include("./physics/energy_emittance_radiative_conductance.jl")
include("./physics/evaporation_vapor_pressure.jl")

# the file modules
include("./file/read_spac_environment.jl")

# the spac modules --- constant
include("./spac/constant.jl")
include("./spac/constant_edit.jl")
include("./spac/constant_edit_soil.jl")
include("./spac/constant_load.jl")
include("./spac/constant_print.jl")

# the spac modules --- energy budget
include("./spac/energy_budget_leaf_temperature.jl")

# the spac modules --- hydraulic supply curve
include("./spac/hydraulic_supply_curve_legacy.jl")
include("./spac/hydraulic_supply_curve_virgin.jl")
include("./spac/hydraulic_supply_etop.jl")
include("./spac/hydraulic_supply_etop_plant.jl")
include("./spac/hydraulic_supply_legacy_initialize.jl")
include("./spac/hydraulic_supply_legacy_print.jl")
include("./spac/hydraulic_supply_legacy_update.jl")
include("./spac/hydraulic_supply_ptoe.jl")

# the spac modules --- photosynthesis
include("./spac/photosynthetic_constants.jl")
include("./spac/photosynthetic_gamma.jl")
include("./spac/photosynthetic_j.jl")
include("./spac/photosynthetic_jmax.jl")
include("./spac/photosynthetic_kc.jl")
include("./spac/photosynthetic_ko.jl")
include("./spac/photosynthetic_photosynthesis.jl")
include("./spac/photosynthetic_vcmax.jl")

# the spac modules --- respiration
include("./spac/respiratic_respiration.jl")

# the spac modules --- simulation
include("./spac/simulation_trade_off.jl")

# the spac modules --- soil
include("./spac/soil_etop.jl")
include("./spac/soil_water_pressure.jl")
include("./spac/soil_water_content.jl")

# the spac modules --- trade-off
include("./spac/trade_off_curves_legacy.jl")
include("./spac/trade_off_curves_virgin.jl")
include("./spac/trade_off_hydraulic_cost.jl")
include("./spac/trade_off_photosynthetic_gain_gross.jl")
include("./spac/trade_off_photosynthetic_gain_net.jl")
include("./spac/trade_off_profit.jl")
include("./spac/trade_off_profit_fix.jl")
include("./spac/trade_off_profit_fix_know_tleaf.jl")
include("./spac/trade_off_profit_optimize.jl")
include("./spac/trade_off_profit_optimize_fix.jl")
include("./spac/trade_off_profit_optimize_fix_know_tleaf.jl")

# the stdio modules
include("./stdio/input_number.jl")
