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
    fmodel,glue=refine(fmodel)
    octree_distributed_discrete_model_free(model)
  end
  if !MPI.Initialized()
    MPI.Init()
  end
  parts = get_part_ids(mpi,1)
  run(parts,(1,1))
  MPI.Finalize()
end # module
