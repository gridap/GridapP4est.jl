module OctreeDistributedDiscreteModelsTests
  using MPI
  using Test
  using Gridap
  using Gridap.ReferenceFEs
  using Gridap.FESpaces
  using PartitionedArrays
  using GridapDistributed
  using GridapP4est
  using P4est_wrapper

  function run(parts,subdomains,num_parts_x_level)
    if length(subdomains) == 2
      domain=(0,1,0,1)
    else
      @assert length(subdomains) == 3
      domain=(0,1,0,1,0,1)
    end

    # Generate model
    level_parts  = GridapP4est.generate_level_parts(parts,num_parts_x_level)
    coarse_model = CartesianDiscreteModel(domain,subdomains)
    model        = OctreeDistributedDiscreteModel(level_parts[2],coarse_model,1)
    vmodel       = GridapP4est._create_void_octree_model(model,parts)

    # Refining and distributing
    fmodel , rglue  = refine(model,level_parts[1])
    dfmodel, dglue  = redistribute(fmodel)

    # FESpaces tests
    sol(x) = x[1] + x[2]
    reffe = ReferenceFE(lagrangian,Float64,1)
    test  = TestFESpace(dfmodel, reffe; conformity=:H1)
    trial = TrialFESpace(sol,test)
 
    # Refining and distributing
    fmodel_tasks_L2, rglue  = refine(model)
    fmodel_tasks_L1, dglueL2toL1  = redistribute(fmodel_tasks_L2,level_parts[1])

    # # Coarsening 
    f_model_tasks_L2_back, dglueL1toL2 = redistribute(fmodel_tasks_L1,level_parts[2])
    # model_back = coarsen(f_model_tasks_L2_back)

    # Destroy models
    octree_distributed_discrete_model_free!(model)
    octree_distributed_discrete_model_free!(fmodel_tasks_L2)
    octree_distributed_discrete_model_free!(fmodel_tasks_L1)
  end

  prun(run,mpi,4,(2,2),[4,2])
  MPI.Finalize()
end # module
