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
    level_parts = Vector{AbstractVector}(undef,nlevs)
    for l = 1:nlevs
      lsize = num_procs_x_level[l]
      if (l > 1) && (lsize == num_procs_x_level[l-1])
        level_parts[l] = level_parts[l-1]
      else
        if rank < lsize
          comm=MPI.Comm_split(root_comm,0,0)
          level_parts[l] = distribute_with_mpi(LinearIndices((lsize,));comm=comm,duplicate_comm=false)
        else
          comm=MPI.Comm_split(root_comm,mpi_undefined,mpi_undefined)
          level_parts[l] = GridapDistributed.MPIVoidVector(eltype(parts))
        end
      end
    end
    return level_parts
  end

  function run_with_env(distribute,np,nc,num_parts_x_level)
    ranks = distribute(LinearIndices((prod(np),)))
    GridapP4est.with(ranks) do
      run(distribute,np,nc,num_parts_x_level)
    end
  end

  function run(distribute,np,nc,num_parts_x_level)
    ranks = distribute(LinearIndices((prod(np),)))
    if length(nc) == 2
      domain=(0,1,0,1)
    else
      @assert length(nc) == 3
      domain=(0,1,0,1,0,1)
    end

    # Generate model
    level_parts  = generate_level_parts(ranks,num_parts_x_level)
    coarse_model = CartesianDiscreteModel(domain,nc)
    model        = OctreeDistributedDiscreteModel(level_parts[2],coarse_model,1) 
    vmodel1      = GridapP4est.VoidOctreeDistributedDiscreteModel(model,
                                                                  ranks)
    vmodel2      = GridapP4est.VoidOctreeDistributedDiscreteModel(coarse_model,
                                                                  ranks,
                                                                  GridapP4est.PXestUniformRefinementRuleType())
    # Refining and distributing
    fmodel , rglue  = refine(model,parts=level_parts[1])
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

    # Redistribute L1 -> L2
    f_model_tasks_L2_back, dglueL1toL2 = redistribute(fmodel_tasks_L1,level_parts[2])

    # Coarsening
    model_back, glue = Gridap.Adaptivity.coarsen(f_model_tasks_L2_back)

    if GridapDistributed.i_am_in(level_parts[2])
      @test num_cells(model_back)==num_cells(model)
      map(model.dmodel.models,model_back.dmodel.models) do m1, m2
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
      omodel, glue = Gridap.Adaptivity.coarsen(imodel)
      imodel = omodel
    end
    @test num_cells(imodel) == prod(nc)
    nothing
  end

  function run_num_ghost_layers_with_env(distribute)
    ranks = distribute(LinearIndices((4,)))
    GridapP4est.with(ranks) do
      test_num_ghost_layers_preserved(ranks)
    end
  end

  # Checks that num_ghost_layers is preserved across refine, adapt, coarsen
  # and redistribute (and across repeated applications of them).
  # A 4x1 strip of coarse cells on 4 ranks, s.t. rank p owns tree p.
  # With each tree refined into n x n cells and g <= n ghost layers,
  # each face-neighbouring tree contributes g*n ghost cells.
  # Ranks 1 and 4 have one neighbouring tree, ranks 2 and 3 have two.
  function test_num_ghost_layers_preserved(ranks)
    @assert length(ranks) == 4
    coarse_model = CartesianDiscreteModel((0,4,0,1),(4,1))
    g = 2

    function test_ghosts(label,model,n)
      @testset "$label" begin
        @test model.num_ghost_layers == g
        map(ranks,partition(get_cell_gids(model))) do rank, indices
          @test ghost_length(indices) == (rank in (1,4) ? 1 : 2)*g*n
        end
      end
    end

    nothing_flags(model) = map(partition(get_cell_gids(model))) do indices
      fill(nothing_flag,local_length(indices))
    end

    @testset "num_ghost_layers preserved" begin
      m0 = OctreeDistributedDiscreteModel(ranks,coarse_model,1;num_ghost_layers=g)
      test_ghosts("m0",m0,2)

      # adapt with all nothing flags (ghost cells do not change)
      a1, _ = Gridap.Adaptivity.adapt(m0,nothing_flags(m0)); test_ghosts("adapt(m0)",a1,2)
      a2, _ = Gridap.Adaptivity.adapt(a1,nothing_flags(a1)); test_ghosts("adapt(adapt(m0))",a2,2)

      # uniform refinement and coarsening (n doubles/halves at each step)
      r1, _ = refine(m0); test_ghosts("refine(m0)",r1,4)
      r2, _ = refine(r1); test_ghosts("refine(refine(m0))",r2,8)
      c1, _ = Gridap.Adaptivity.coarsen(r2); test_ghosts("coarsen(r2)",c1,4)
      c2, _ = Gridap.Adaptivity.coarsen(c1); test_ghosts("coarsen(coarsen(r2))",c2,2)

      # redistribution of an already balanced mesh (ghost cells do not change)
      d1, _ = redistribute(m0); test_ghosts("redistribute(m0)",d1,2)
      d2, _ = redistribute(d1); test_ghosts("redistribute(redistribute(m0))",d2,2)
    end

    # Void models (on ranks not in the model's parts) must also carry
    # num_ghost_layers, as they take part in collective ghost layer
    # construction when refining/redistributing onto a larger set of parts.
    @testset "num_ghost_layers preserved by void models" begin
      level_parts = generate_level_parts(ranks,[4,2])
      vm1 = GridapP4est.VoidOctreeDistributedDiscreteModel(coarse_model,ranks,
              GridapP4est.PXestUniformRefinementRuleType();num_ghost_layers=g)
      @test vm1.num_ghost_layers == g

      # Void on ranks 3 and 4
      m2 = OctreeDistributedDiscreteModel(level_parts[2],coarse_model,1;num_ghost_layers=g)
      @test m2.num_ghost_layers == g
      vm2 = GridapP4est.VoidOctreeDistributedDiscreteModel(m2,ranks)
      @test vm2.num_ghost_layers == g

      # 2 -> 4 parts: rank p owns tree p
      d4, _ = redistribute(m2,ranks); test_ghosts("redistribute(m2,ranks)",d4,2)
      # 4 -> 2 parts: void on ranks 3 and 4
      d2, _ = redistribute(d4,level_parts[2])
      @test d2.num_ghost_layers == g
      # 2 -> 4 parts again, starting from the void models above
      d4b, _ = redistribute(d2,ranks); test_ghosts("redistribute(d2,ranks)",d4b,2)

      # refine onto 4 parts, then balance s.t. rank p owns tree p
      r4, _ = refine(m2;parts=ranks)
      @test r4.num_ghost_layers == g
      r4b, _ = redistribute(r4); test_ghosts("redistribute(refine(m2;parts=ranks))",r4b,4)
    end
    nothing
  end
end # module
