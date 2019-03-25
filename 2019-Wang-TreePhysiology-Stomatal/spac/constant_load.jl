# This function is used to load known experimental data so as to save time, currently the data loading are for the water birch trees.

# BOS stand for Betula occidentalis in 2017
function ConstantLoadBO2017()
    tmp_parts = Array[
        [367.3476, 1.56, 0.358974, 3E13],
        [1.8794, 2.3956,  4657.32],
        [2.2377, 9.3796, 10298.08],
        [1.8968, 2.2026, 11338.24 ,107.814, 1.28196, 50.0, 400.0, 20000.0, 0.05, 0.5, 10.0, 0.001, 4E5, 5E-4]]
    return tmp_parts
end
