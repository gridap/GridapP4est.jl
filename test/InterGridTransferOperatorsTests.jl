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

  function change_domain_coarse_to_fine(c_cell_field,
                                        ftrian::Triangulation{Dc,Dp},
                                        glue::FineToCoarseModelGlue) where {Dc,Dp}

    fcell_to_child_id=glue.fcell_to_child_id
    fcell_to_ccell=glue.fine_to_coarse_faces_map[Dc+1]

    Gridap.Helpers.@check length(fcell_to_child_id)==num_cells(ftrian)
    Gridap.Helpers.@check DomainStyle(c_cell_field)==ReferenceDomain()

    rrule_reffe=Gridap.ReferenceFEs.LagrangianRefFE(Float64,QUAD,1)
    rrule_grid=Gridap.Geometry.UnstructuredGrid(
                  Gridap.Visualization.compute_reference_grid(rrule_reffe,2))

    rrule_f_to_c_ref_map=get_cell_map(rrule_grid)
    cfield=Gridap.CellData.get_data(c_cell_field)
    m1=Reindex(cfield)
    f_to_cfield=lazy_map(m1,fcell_to_ccell)
    m2=Reindex(rrule_f_to_c_ref_map)
    cfield_to_ffield=lazy_map(∘,f_to_cfield,lazy_map(m2,fcell_to_child_id))
    Gridap.CellData.GenericCellField(cfield_to_ffield,ftrian,ReferenceDomain())
  end

  function change_domain_coarse_to_fine(c_cell_field::GridapDistributed.DistributedCellField,
                                        ftrian::GridapDistributed.DistributedTriangulation{Dc,Dp},
                                        glue::MPIData{<:FineToCoarseModelGlue}) where {Dc,Dp}

    dfield=map_parts(change_domain_coarse_to_fine,
                     GridapDistributed.local_views(c_cell_field),
                     GridapDistributed.local_views(ftrian),
                     glue)
    GridapDistributed.DistributedCellField(dfield)
  end



  function ProlongationOperator(uH,Uh,Vh,ftrian,glue,dΩ)
    uH_h = change_domain_coarse_to_fine(uH,ftrian,glue)
    a(u,v) = ∫(v*u)dΩ
    l(v)   = ∫(v*uH_h)dΩ
    op=AffineFEOperator(a,l,Uh,Vh)
    solve(op)
  end


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
    uH_h   = change_domain_coarse_to_fine(uH,ftrian,glue)

    trian=Triangulation(fmodel.dmodel)
    dΩ=Measure(trian,2*(order+1))

    uh=ProlongationOperator(uH,Uh,Vh,trian,glue,dΩ)

    e = uH_h-uh
    e_l2 = sum(∫(e*e)dΩ)
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
