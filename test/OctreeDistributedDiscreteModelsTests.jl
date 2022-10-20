module OctreeDistributedDiscreteModelsTests
  using MPI
  using Gridap
  using Gridap.ReferenceFEs
  using Gridap.FESpaces
  using PartitionedArrays
  using GridapDistributed
  using GridapP4est
  using P4est_wrapper

  function run(parts,subdomains)
    if length(subdomains) == 2
      domain=(0,1,0,1)
    else
      @assert length(subdomains) == 3
      domain=(0,1,0,1,0,1)
    end

    # Generate models
    coarse_model = CartesianDiscreteModel(domain,subdomains)
    model        = OctreeDistributedDiscreteModel(parts,coarse_model,1)
    fmodel,glue  = refine(model)
    ffmodel,glue = refine(fmodel)

    # FESpaces tests
    sol(x) = x[1] + x[2]
    reffe = ReferenceFE(lagrangian,Float64,1)
    test  = TestFESpace(model, reffe; conformity=:H1)
    trial = TrialFESpace(sol,test)

    # Destroy models
    octree_distributed_discrete_model_free!(model)
    octree_distributed_discrete_model_free!(fmodel)
    octree_distributed_discrete_model_free!(ffmodel)
  end

  prun(run,mpi,4,(2,2))
  MPI.Finalize()
end # module
