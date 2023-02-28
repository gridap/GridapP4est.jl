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

  import Gridap.Adaptivity: refine
  import GridapDistributed: redistribute

  function generate_level_parts(parts,num_procs_x_level)
    root_comm = parts.comm
    rank = MPI.Comm_rank(root_comm)
    size = MPI.Comm_size(root_comm)
    Gridap.Helpers.@check all(num_procs_x_level .<= size)
    Gridap.Helpers.@check all(num_procs_x_level .>= 1)
  
    @static if isdefined(MPI,:MPI_UNDEFINED)
      mpi_undefined = MPI.MPI_UNDEFINED[]
    else
      mpi_undefined = MPI.API.MPI_UNDEFINED[]
    end
    
    nlevs = length(num_procs_x_level)
    level_parts = Vector{typeof(parts)}(undef,nlevs)
    for l = 1:nlevs
      lsize = num_procs_x_level[l]
      if (l > 1) && (lsize == num_procs_x_level[l-1])
        level_parts[l] = level_parts[l-1]
      else
        if rank < lsize
          comm=MPI.Comm_split(root_comm,0,0)
        else
          comm=MPI.Comm_split(root_comm,mpi_undefined,mpi_undefined)
        end
        level_parts[l] = get_part_ids(comm)
      end
    end
    return level_parts
  end

  function run_with_env(parts,subdomains,num_parts_x_level)
    GridapP4est.with(parts) do
      run(parts,subdomains,num_parts_x_level)
    end
  end

  function run(parts,subdomains,num_parts_x_level)
    if length(subdomains) == 2
      domain=(0,1,0,1)
    else
      @assert length(subdomains) == 3
      domain=(0,1,0,1,0,1)
    end

    # Generate model
    level_parts  = generate_level_parts(parts,num_parts_x_level)
    coarse_model = CartesianDiscreteModel(domain,subdomains)
    model        = OctreeDistributedDiscreteModel(level_parts[2],coarse_model,1)
    vmodel1      = GridapP4est.VoidOctreeDistributedDiscreteModel(model,parts)
    vmodel2      = GridapP4est.VoidOctreeDistributedDiscreteModel(coarse_model,parts)

    # Refining and distributing
    fmodel , rglue  = refine(model,level_parts[1])
    dfmodel, dglue  = redistribute(fmodel)

    # FESpaces tests
    sol(x) = x[1] + x[2]
    reffe = ReferenceFE(lagrangian,Float64,1)
    test  = TestFESpace(dfmodel, reffe; conformity=:H1)
    trial = TrialFESpace(sol,test)

    # Refine
    fmodel_tasks_L2, rglue  = refine(model)

    # Redistribute L2 -> L1
    fmodel_tasks_L1, dglueL2toL1  = redistribute(fmodel_tasks_L2,level_parts[1])
    if i_am_in(level_parts[1])
      @test fmodel_tasks_L1.parts === PartitionedArrays.get_part_ids(dglueL2toL1.parts_rcv)
    end

    # Redistribute L1 -> L2
    f_model_tasks_L2_back, dglueL1toL2 = redistribute(fmodel_tasks_L1,level_parts[2])

    # Coarsening
    model_back,glue = coarsen(f_model_tasks_L2_back)

    if i_am_in(level_parts[2])
      @test num_cells(model_back)==num_cells(model)
      map_parts(model.dmodel.models,model_back.dmodel.models) do m1, m2
        Ωh1  = Triangulation(m1)
        dΩh1 = Measure(Ωh1,2)
        Ωh2  = Triangulation(m2)
        dΩh2 = Measure(Ωh2,2)
        sum(∫(1)dΩh1) ≈ sum(∫(1)dΩh2)
      end
    end

    model  = OctreeDistributedDiscreteModel(level_parts[1],coarse_model,3)
    imodel = model
    for i=1:3
      omodel,glue=coarsen(imodel)
      imodel=omodel
    end
    @test num_cells(imodel)==prod(subdomains)
    nothing
  end
end # module
