# Published Codes by Yujie WANG

## Name Rules
To make it easy for readers to locate the models, the projects are named by the rule of Year-Author-Journal-FirstTitleWord.

## Projects
| Folder                                | Description                                                                               |
|:--------------------------------------|:------------------------------------------------------------------------------------------|
| 2024-Wang-GCB-Beyond                  | Julia and Python code used for data analysis and plotting                                 |
| 2024-Wang-JAMES-Toward                | Julia and Python code used for data analysis and plotting                                 |
| 2023-Wang-CROPE-Agriculture           | Julia code for the prototype models for crop studies                                      |
| 2023-Braghiere-AGUAdvances-Importance | Julia code for the hyperspectral soil albedo project                                      |
| 2023-Wang-JAMES-Modeling              | Julia code used to plot the figures used in CliMA Land global simulations                 |
| 2022-Wang-BG-Common                   | Julia code used to plot the figure for plant hydraulics models                            |
| 2022-Wang-SDATA-GriddingMachine       | Julia code used to plot the figure for GriddingMachine illustration                       |
| 2022-Wang-BG-Impact                   | Julia code used to evluate the impact of canopy model complexity on water and SIF fluxes  |
| 2021-Wang-GMD-Testing                 | Julia code used to run site level CliMA Land simulations                                  |
| 2021-Braghiere-RSE-Accounting         | Julia code used to fit leaf biophysical parameters                                        |
| 2021-Wang-NPh-Optimization            | Julia code used to simulate optimal nighttime stomatal conductance                        |
| 2020-Wang-NPh-Theoretical             | Python code for stomatal optimization models                                              |
| 2019-Wang-TPh-Stomatal                | Julia version of Sperry gain-risk model (basic version)                                   |
| 2018-Venturas-NPh-Stomatal            | Excel version of Sperry gain-risk model (more updated)                                    |
| 2017-Sperry-PCE-Predicting            | Excel version of Sperry gain-risk model                                                   |
| 2016-Sperry-NPh-Pragmatic             | C and Python cide to combine the VCs of root, stem, and leaf into an integrated tree VC   |
| 2015-Wang-PPh-Stem                    | Python code to fit bubble pressure in conduits of recently cavitated stem from two points |
| 2015-Wang-PPh-Studies                 | Python code to fit bubble pressure in conduits of recently cavitated stem from a curve    |
| 2014-Wang-JPH-Improved                | Software for Chinatron and brief manual                                                   |
|||

## Note
For the Julia projects, the Manifest.toml files were provided for the exact model version as the models are actively being developed. Therefore, you need to use
```julia
pkg> instantiate
```
Otherwise, the code likely would not work if you update to the latest code.
