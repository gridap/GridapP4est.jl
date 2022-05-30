module PCMGTests
  using MPI
  using Gridap
  using GridapPETSc
  using PartitionedArrays
  using GridapDistributed
  using GridapP4est
  using P4est_wrapper
  using Test

  # Manufactured solution
  u(x) = x[1] + x[2]
  f(x) = -Δ(u)(x)

  function run(parts,subdomains)
    # if length(subdomains)==2
    #   domain=(0,1,0,1)
    # else
    #   @assert length(subdomains)==3
    #   domain=(0,1,0,1,0,1)
    # end
    # coarse_discrete_model=CartesianDiscreteModel(domain,subdomains)
    # model=OctreeDistributedDiscreteModel(parts,
    #                                      coarse_discrete_model)
    # cmodel,_=refine(model)
    # fmodel,glue=refine(cmodel)

    # # FE Spaces
    # order=1
    # reffe=ReferenceFE(lagrangian,Float64,order)
    # VH=TestFESpace(cmodel.dmodel,reffe,dirichlet_tags="boundary")
    # UH=TrialFESpace(u,VH)
    # Vh=TestFESpace(fmodel.dmodel,reffe,dirichlet_tags="boundary")
    # Uh=TrialFESpace(u,Vh)

    # uH     = interpolate(u,UH)
    # ftrian = get_triangulation(fmodel.dmodel)
    # ctrian = get_triangulation(cmodel.dmodel)
    # uH_h   = GridapP4est.change_domain_coarse_to_fine(uH,ftrian,glue)

    # dΩh=Measure(ftrian,2*(order+1))
    # dΩH=Measure(ctrian,2*(order+1))

    # uh=ProlongationOperator(uH,Uh,Vh,ftrian,glue,dΩh)
    # uH_from_uh=RestrictionOperator(uh,UH,VH,ctrian,glue,dΩH)


    function setup_KSP_PCMG(comm::MPI.Comm)
      ksp=Ref{KSP}()
      pc=Ref{PC}()
      @check_error_code GridapPETSc.PETSC.KSPCreate(comm,ksp)
      @check_error_code GridapPETSc.PETSC.KSPSetFromOptions(ksp[])
      @check_error_code GridapPETSc.PETSC.KSPGetPC(ksp[],pc)
      @check_error_code GridapPETSc.PETSC.PCSetType(pc[],GridapPETSc.PETSC.PCMG)
      @check_error_code GridapP4est.PCMGSetLevels(pc[],PetscInt(2),C_NULL)
      @check_error_code GridapPETSc.PETSC.KSPView(ksp[],C_NULL)
      ksp
    end

    options = "-ksp_type cg -pc_type mg -pc_mg_type multiplicative -pc_mg_multiplicative_cycles 1 -mg_coarse_pc_type cholesky -mg_coarse_ksp_type preonly -mg_levels_pc_type cholesky -mg_levels_ksp_type preonly"
    GridapPETSc.with(args=split(options)) do
       ksp=setup_KSP_PCMG(MPI.COMM_WORLD)
       GridapPETSc.PETSC.KSPDestroy(ksp)
    end
    # octree_distributed_discrete_model_free(cmodel)
    # octree_distributed_discrete_model_free(fmodel)
  end
  if !MPI.Initialized()
    MPI.Init()
  end
  parts = get_part_ids(mpi,1)
  run(parts,(1,1))
  MPI.Finalize()
end
