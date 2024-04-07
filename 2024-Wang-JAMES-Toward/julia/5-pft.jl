#
# this script is meant for compute the doninant PFT of each grid cell
#
using GriddingMachine.Blender: regrid
using GriddingMachine.Collector: query_collection
using GriddingMachine.Indexer: read_LUT
using NetcdfIO: save_nc!


# 1. read the PFT data
pft_map = regrid(read_LUT(query_collection("PFT_2X_1Y_V1"))[1], 1);
CLM5_PFTS = ["not_vegetated",
             "needleleaf_evergreen_temperate",
             "needleleaf_evergreen_boreal",
             "needleleaf_deciduous_boreal",
             "broadleaf_evergreen_tropical",
             "broadleaf_evergreen_temperate",
             "broadleaf_deciduous_tropical",
             "broadleaf_deciduous_temperate",
             "broadleaf_deciduous_boreal",
             "evergreen_shrub",
             "deciduous_temperate_shrub",
             "deciduous_boreal_shrub",
             "c3_arctic_grass",
             "c3_non-arctic_grass",
             "c4_grass",
             "c3_crop",
             "c3_irrigated"];

# 2. choose only these PFTs: 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12
dominating_pft = zeros(Int, size(pft_map,1), size(pft_map,2));
pft_choices = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];
pft_counts = zeros(Int, length(pft_choices));
for ilon in eachindex(axes(pft_map, 1))
    for ilat in eachindex(axes(pft_map, 2))
        pft_fractions = pft_map[ilon, ilat, :][pft_choices];
        if sum(pft_fractions) > 0
            max_i = findmax(pft_fractions)[2];
            pft_counts[max_i] += 1;
            max_pft_ind = pft_choices[max_i];
            dominating_pft[ilon, ilat] = max_pft_ind;
        end;
    end;
end;

# 3. save the data
save_nc!("../output/5_pft.nc", "PFT", dominating_pft, Dict{String,String}("about" => "Dominating PFT in each grid cell"));
