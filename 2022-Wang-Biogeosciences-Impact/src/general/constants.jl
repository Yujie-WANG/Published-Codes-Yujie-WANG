###############################################################################
#
# Constant strings
#
###############################################################################
const LS_Anet   = latex_symbol("A", sub="net"  );
const LS_gsw    = latex_symbol("g", sub="sw"   );
const LS_Vcmax  = latex_symbol("V", sub="cmax" );

const LS_unit_SIF     = "(mW m\$^{-2}\$ sr\$^{-1}\$ nm\$^{-1}\$)";

const LS_APAR_unit   = "APAR " * latex_unit("PAR");
const LS_Anet_unit   = LS_Anet * " " * latex_unit("A");
const LS_ET_unit_m   = "ET " * latex_unit("E_MMOL");
const LS_gsw_unit    = LS_gsw * " " * latex_unit("G");
const LS_NEE_unit    = "NEE " * latex_unit("A");
const LS_SIF740_unit = "SIF\$_{740}\$ $(LS_unit_SIF)";
const LS_Vcmax_unit  = LS_Vcmax * " " * latex_unit("A");


###############################################################################
#
# Constant dictionaries
#
###############################################################################
const ARROW_C  = Dict(["arrowstyle" => "->" , "color" => "c"]);
const ARROW_G  = Dict(["arrowstyle" => "->" , "color" => "g"]);
const COLORS   = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"];
