###############################################################################
#
# Dynamically change number of workers
#
# Move this to PlotPlants?
#
###############################################################################
"""
Add or remove threads dynamically.
"""
function dynamic_workers! end




"""
This method adds or removes workers to a certain number. Note that is function
    add workers only, but does not load anything into the worker. Given that
    the workers are meant to work under the developing environment, all of them
    have a `--project` flag.

    dynamic_workers!(nTH::Int)

Change the number of workers according to `CPU_THREADS`, given
- `nTH` Number of thread required
"""
dynamic_workers!(nTH::Int) =
(
    _MaxThreads = Sys.CPU_THREADS;

    # if        no worker yet, and nTH <= MaxThreads, hire nTH
    # elseif    no worker yet, but nTH > MaxThreads, hire MaxThreads
    # elseif    some workers already, and nTH <= MaxThreads, hire more
    # elseif    some workers already, and nTH > MaxThreads, hire more
    # else      workers is more than expected, remove the extra
    if (length(workers())==1) && (workers()[1]==1) && (nTH<=_MaxThreads)
        addprocs(nTH, exeflags="--project");
    elseif (length(workers())==1) && (workers()[1]==1) && (nTH>_MaxThreads)
        addprocs(_MaxThreads, exeflags="--project");
    elseif length(workers())<nTH && (nTH<=_MaxThreads)
        addprocs(nTH-length(workers()), exeflags="--project");
    elseif length(workers())<_MaxThreads && (nTH>_MaxThreads)
        addprocs(_MaxThreads-length(workers()), exeflags="--project");
    else
        _to_remove = workers()[(nTH+1):end];
        rmprocs(_to_remove...);
    end;

    return nothing
)




"""
This method builds on the method above and loads required packages to threads:

    dynamic_workers!(
                project::ForestTower2021{FT},
                nTH::Int
    ) where {FT<:AbstractFloat}

Change the number of workers according to `CPU_THREADS`, given
- `proj` Project identifier. Is one of the following
  - [`ForestTower2021`](@ref)
  - [`LandGPP2021`](@ref)
  - [`SIFGPP2021`](@ref)
- `nTH` Number of thread required
"""
dynamic_workers!(
            project::ForestTower2021{FT},
            nTH::Int
) where {FT<:AbstractFloat} =
(
    dynamic_workers!(nTH);

    # load the module into each worker
    @everywhere Base.MainInclude.eval(using ResearchProjects);

    return nothing
)
