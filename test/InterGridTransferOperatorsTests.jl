module InterGridTransferOperatorsTests
  using MPI
  using Gridap
  using PartitionedArrays
  using GridapDistributed
  using GridapP4est
  using P4est_wrapper
  using Test

  # Manufactured solution
  u(x) = x[1] + x[2]
  f(x) = -Δ(u)(x)

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
    cmodel,_=refine(model)
    fmodel,glue=refine(cmodel)

    # FE Spaces
    order=1
    reffe=ReferenceFE(lagrangian,Float64,order)
    VH=TestFESpace(cmodel.dmodel,reffe,dirichlet_tags="boundary")
    UH=TrialFESpace(u,VH)
    Vh=TestFESpace(fmodel.dmodel,reffe,dirichlet_tags="boundary")
    Uh=TrialFESpace(u,Vh)

    uH     = interpolate(u,UH)
    ftrian = get_triangulation(fmodel.dmodel)
    ctrian = get_triangulation(cmodel.dmodel)
    uH_h   = GridapP4est.change_domain_coarse_to_fine(uH,ftrian,glue)

    dΩh=Measure(ftrian,2*(order+1))
    dΩH=Measure(ctrian,2*(order+1))

    uh=ProlongationOperator(uH,Uh,Vh,ftrian,glue,dΩh)
    uH_from_uh=RestrictionOperator(uh,UH,VH,ctrian,glue,dΩH)

    e = uH_h-uh
    e_l2 = sum(∫(e*e)dΩh)
    tol = 1.0e-9
    @test e_l2 < tol
    map_parts(parts) do part
     if (part==1)
       println("$(e_l2) < $(tol)\n")
     end
    end

    e = uH_from_uh-uH
    e_l2 = sum(∫(e*e)dΩH)
    tol = 1.0e-9
    @test e_l2 < tol
    map_parts(parts) do part
     if (part==1)
       println("$(e_l2) < $(tol)\n")
     end
    end

    octree_distributed_discrete_model_free(cmodel)
    octree_distributed_discrete_model_free(fmodel)
  end
  if !MPI.Initialized()
    MPI.Init()
  end
  parts = get_part_ids(mpi,1)
  run(parts,(1,1))
  MPI.Finalize()
end
