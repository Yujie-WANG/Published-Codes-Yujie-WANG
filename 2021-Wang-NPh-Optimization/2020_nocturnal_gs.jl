using ResearchProjects




FT     = Float32;
proj   = NocturnalGS2020{FT}();
saving = true;

@info "Now plotting Fig. 1...";
plot_model_framework(proj; saving=saving);

@info "Now plotting Fig. 2...";
plot_model_prediction(proj; saving=saving);

@info "Now plotting Fig. 3-5...";
plot_model_comparison(proj; saving=saving);

@info "Now plotting Fig. 6...";
plot_gswn_vs_time(proj; saving=saving);

@info "Now plotting Fig. 7...";
plot_model_extension(proj; saving=saving);

@info "Now plotting Fig. S1...";
plot_si_t_leaf(proj; saving=saving);
