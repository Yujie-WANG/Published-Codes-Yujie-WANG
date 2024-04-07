### A Pluto.jl notebook ###
# v0.17.5

using Markdown
using InteractiveUtils

# ╔═╡ 6b00d686-26f0-11ec-2140-2bb8bb37199b
# load packages
begin
    using Pkg;
    Pkg.activate("Project.toml");

    using PlotPlants: create_canvas, save_canvas!, set_xlabels!, set_xlims!, set_ylabels!, set_ylims!;
end;

# ╔═╡ d8cccc67-b2d3-4303-a974-abd23d68511b
# request data from the server for US-NR1
begin
    using GriddingMachine.Requestor;
    # load FLUXCOM and VPM GPPs and TROPOMI DC SIF
    gpp_mpi,_ = request_LUT("GPP_MPI_RS_2X_8D_2019_V1", 40.0329, -105.5464; user="Yujie");
    gpp_vpm,_ = request_LUT("GPP_VPM_12X_8D_2019_V2", 40.0329, -105.5464; user="Yujie");
    sif_dc,_  = request_LUT("SIF_TROPOMI_740DC_12X_8D_2019_V1", 40.0329, -105.5464; user="Yujie");
end;

# ╔═╡ f2be78a6-f43b-4086-9832-eddc105bdbd2
# set width to 1300 px
html"""<style>
main {
    max-width: 1300px;
    margin-right: 350px;
}
"""

# ╔═╡ d297f803-df14-45fe-a791-0dccf21c3918
# plot the results from the server for US-NR1
begin
    # plot the comparison
    fig,axs = create_canvas("GM-Example for US-NR1", figsize=(6.5,3.5));
    ax1, = axs;
    tx1 = ax1.twinx();
    ax1.plot(4:8:365, gpp_mpi, "yo-", label="MPI GPP");
    ax1.plot(4:8:365, gpp_vpm, "co-", label="VPM GPP");
    tx1.plot(4:8:365, sif_dc, "ro:", mfc="none", label="TROPOMI SIF");
    ax1.legend(loc="upper left");
    tx1.legend(loc="upper right");
    set_xlims!(axs, [0,368]);
    set_ylims!([ax1,tx1], [[-2.2,5.5], [-0.1,0.25]]);
    set_xlabels!(axs, "Day of year 2019");
    set_ylabels!([ax1,tx1], ["GPP (g C m⁻² day⁻¹)", "dcSIF (mW m² nm⁻¹ sr⁻¹)"]);
    save_canvas!(fig, "3_example.pdf", true);
    fig;
end

# ╔═╡ Cell order:
# ╠═f2be78a6-f43b-4086-9832-eddc105bdbd2
# ╠═6b00d686-26f0-11ec-2140-2bb8bb37199b
# ╠═d8cccc67-b2d3-4303-a974-abd23d68511b
# ╠═d297f803-df14-45fe-a791-0dccf21c3918
