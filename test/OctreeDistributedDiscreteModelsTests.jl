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

  function run(parts,subdomains,num_parts_x_level)
    GridapP4est.with(parts) do
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
      if GridapP4est.i_am_in(level_parts[1])
        @test fmodel_tasks_L1.parts === PartitionedArrays.get_part_ids(dglueL2toL1.parts_rcv)
      end

      # Redistribute L1 -> L2
      f_model_tasks_L2_back, dglueL1toL2 = redistribute(fmodel_tasks_L1,level_parts[2])

      # Coarsening
      model_back,glue = coarsen(f_model_tasks_L2_back)

      if (GridapP4est.i_am_in(level_parts[2]))
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
  end

  prun(run,mpi,4,(2,2),[4,2])
  MPI.Finalize()
end # module
