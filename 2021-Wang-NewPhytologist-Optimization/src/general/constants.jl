###############################################################################
#
# Constant strings
#
###############################################################################
const LS_Ad    = latex_symbol("A", sub="d"    );
const LS_An    = latex_symbol("A", sub="n"    );
const LS_Anet  = latex_symbol("A", sub="net"  );
const LS_Ca    = latex_symbol("C", sub="a"    );
const LS_Cc    = latex_symbol("C", sub="c"    );
const LS_Ci    = latex_symbol("C", sub="i"    );
const LS_Ed    = latex_symbol("E", sub="d"    );
const LS_En    = latex_symbol("E", sub="n"    );
const LS_ff    = latex_symbol("f", sub="f"    );
const LS_gm    = latex_symbol("g", sub="m"    );
const LS_gmax  = latex_symbol("g", sub="max"  );
const LS_gs    = latex_symbol("g", sub="s"    );
const LS_gsc   = latex_symbol("g", sub="sc"   );
const LS_gwn   = latex_symbol("g", sub="wn"   );
const LS_gwn25 = latex_symbol("g", sub="wn,25");
const LS_Pleaf = latex_symbol("P", sub="leaf" );
const LS_Kx    = latex_symbol("K", sub="x"    );
const LS_Psoil = latex_symbol("P", sub="soil" );
const LS_Rbase = latex_symbol("R", sub="base" );
const LS_Rleaf = latex_symbol("R", sub="leaf" );
const LS_Tleaf = latex_symbol("T", sub="leaf" );
const LS_Vcmax = latex_symbol("V", sub="cmax" );

const LS_Ψbase      = "\$\\Uppsi_\\text{base}\$";
const LS_∂A∂gm      = "\$\\partial A / \\partial g_\\text{m}\$";
const LS_∂A∂gm_frac = "\$\\dfrac{\\partial A}{\\partial g_\\text{m}}\$";
const LS_ff_∂Ad∂Ed  = "\$f_\\mathrm{f} \\cdot " *
                      "\\partial A_\\mathrm{d} / \\partial E_\\mathrm{d}\$";
const LS_ff_∂Θd∂Ed  = "\$f_\\mathrm{f} \\cdot " *
                      "\\partial \\Theta_\\mathrm{d} /" *
                      "\\partial E_\\mathrm{d}\$";
const LS_∂Rn∂En     = "\$\\partial R_\\mathrm{leaf} / \\partial E_\\mathrm{n}\$";
const LS_n_∂Rn∂En   = "\$-\\partial R_\\mathrm{leaf} / \\partial E_\\mathrm{n}\$";

const LS_Anet_unit  = LS_Anet * " " * latex_unit("A");
const LS_Ca_unit    = LS_Ca * " (ppm)";
const LS_Cc_unit    = LS_Cc * " (ppm)";
const LS_Ci_unit    = LS_Ci * " (ppm)";
const LS_E_unit     = "\$E\$ " * latex_unit("E");
const LS_E_unit_m   = "\$E\$ " * latex_unit("E_MMOL");
const LS_gm_unit    = LS_gm * " " * latex_unit("G");
const LS_gmax_unit  = LS_gmax * " " * latex_unit("G");
const LS_GPP_unit   = "GPP " * latex_unit("A");
const LS_gwn_unit   = LS_gwn * " " * latex_unit("G");
const LS_gwn25_unit = LS_gwn25 * " " * latex_unit("G");
const LS_kx_unit    = LS_Kx * " (mol H\$_2\$O s\$^{-1}\$ MPa\$^{-1}\$)";
const LS_NEE_unit   = "NEE " * latex_unit("A");
const LS_PAR_unit   = "PAR " * latex_unit("PAR");
const LS_Pleaf_unit = LS_Pleaf * " (MPa)";
const LS_Psoil_unit = LS_Psoil * " (MPa)";
const LS_Rleaf_unit = LS_Rleaf * " " * latex_unit("A");
const LS_Rbase_unit = LS_Rbase * " " * latex_unit("A");
const LS_SIF_unit   = "SIF (mW m\$^{-2}\$ sr\$^{-1}\$ nm\$^{-1}\$)";
const LS_Tleaf_unit = LS_Tleaf * " (\$^{\\circ}\$C)";
const LS_Vcmax_unit = LS_Vcmax * " " * latex_unit("A");
const LS_Ψbase_unit = LS_Ψbase * " (MPa)";
const LS_∂A∂gm_unit = LS_∂A∂gm * latex_unit("WUE");








###############################################################################
#
# Constant dictionaries
#
###############################################################################
const ARROW_B  = Dict(["arrowstyle" => "->" , "color" => "b"]);
const ARROW_C  = Dict(["arrowstyle" => "->" , "color" => "c"]);
const ARROW_G  = Dict(["arrowstyle" => "->" , "color" => "g"]);
const ARROW_K  = Dict(["arrowstyle" => "->" , "color" => "k"]);
const ARROW_R  = Dict(["arrowstyle" => "->" , "color" => "r"]);
const DARROW_K = Dict(["arrowstyle" => "<->", "color" => "k"]);
