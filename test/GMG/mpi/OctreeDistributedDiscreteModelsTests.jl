module OctreeDistributedDiscreteModelsTests
  using MPI
  using Gridap
  using PartitionedArrays
  using GridapDistributed
  using GridapP4est
  using P4est_wrapper

  function run(parts,subdomains)
    if length(subdomains)==2
      domain=(0,1,0,1)
    else
      @assert length(subdomains)==3
      domain=(0,1,0,1,0,1)
    end
    coarse_discrete_model=CartesianDiscreteModel(domain,subdomains)
    model=OctreeDistributedDiscreteModel(parts,
                                         coarse_discrete_model)
    fmodel,glue=refine(model)
    ffmodel,glue=refine(fmodel)
    octree_distributed_discrete_model_free!(model)
    octree_distributed_discrete_model_free!(fmodel)
    octree_distributed_discrete_model_free!(ffmodel)
  end
  if !MPI.Initialized()
    MPI.Init()
  end
  parts = get_part_ids(mpi,6)
  run(parts,(2,2))
  MPI.Finalize()
end # module
