module RedistributeToolsTests
  using MPI
  using PartitionedArrays
  using Gridap
  using GridapDistributed
  using GridapP4est
  using P4est_wrapper
  using Test

  function Base.map(::typeof(Gridap.Arrays.testitem),
                  a::Tuple{<:AbstractVector{<:AbstractVector{<:VectorValue}},<:AbstractVector{<:Gridap.Fields.LinearCombinationFieldVector}})
    a2=Gridap.Arrays.testitem(a[2])
    a1=Vector{eltype(eltype(a[1]))}(undef,size(a2,1))
    a1.=zero(Gridap.Arrays.testitem(a1))
    (a1,a2)
  end

  # Fixes Err3 (see below)
  function Gridap.Geometry.is_change_possible(
      strian::Gridap.Geometry.Triangulation,
      ttrian::Gridap.Geometry.Triangulation)
    if strian === ttrian || num_cells(strian)==num_cells(ttrian)==0
      return true
    end
    Gridap.Helpers.@check get_background_model(strian) === get_background_model(ttrian) "Triangulations do not point to the same background discrete model!"
    D = num_cell_dims(strian)
    sglue = get_glue(strian,Val(D))
    tglue = get_glue(ttrian,Val(D))
    Gridap.Geometry.is_change_possible(sglue,tglue) # Fails here
  end

  # Fixes Err3 (see below)
  function Gridap.CellData.change_domain(a::CellField,
                                        ::ReferenceDomain,
                                        ttrian::Triangulation,
                                        ::ReferenceDomain)
    msg = """\n
    We cannot move the given CellField to the reference domain of the requested triangulation.
    Make sure that the given triangulation is either the same as the triangulation on which the
    CellField is defined, or that the latter triangulation is the background of the former.
    """
    strian = get_triangulation(a)
    if strian === ttrian || num_cells(strian)==num_cells(ttrian)==0
      return a
    end
    @assert Gridap.Geometry.is_change_possible(strian,ttrian) msg
    D = num_cell_dims(strian)
    sglue = get_glue(strian,Val(D))
    tglue = get_glue(ttrian,Val(D))
    Gridap.CellData.change_domain_ref_ref(a,ttrian,sglue,tglue)
  end

  u(x) = x[1] + x[2]

  function run(parts,num_parts_x_level)
    domain=(0,1,0,1)
    partition=(1,1)
    cmodel=CartesianDiscreteModel(domain,partition)
    mh=ModelHierarchy(parts,cmodel,num_parts_x_level)

    # FE Spaces
    order=1
    reffe = ReferenceFE(lagrangian,Float64,order)

    VOLD  = TestFESpace(mh.levels[1].model.dmodel,reffe,dirichlet_tags="boundary")

    # map_parts(mh.level_parts[1],mh.levels[1].model.dmodel.models) do part, model
    #   if (part==2)
    #     Gridap.Io.to_json_file(model,"model_part_2")
    #     println(model)
    #   end
    # end

    UOLD  = TrialFESpace(u,VOLD)
    uhold = interpolate(u,UOLD)

    VNEW = TestFESpace(mh.levels[1].model_red.dmodel,reffe,dirichlet_tags="boundary")
    UNEW = TrialFESpace(u,VNEW)

    uhnew  = redistribute_fe_function(uhold,
                                      UOLD,
                                      UNEW,
                                      mh.levels[1].model_red.dmodel,
                                      mh.levels[1].red_glue)

    Ω_old  = Triangulation(mh.levels[1].model.dmodel)
    dΩ_old = Measure(Ω_old,2)

    Ω_new  = Triangulation(mh.levels[1].model_red.dmodel)
    dΩ_new = Measure(Ω_new,2)

    o=sum(∫(uhold)dΩ_old)
    n=sum(∫(uhnew)dΩ_new)

    @test o ≈ n

    # map_parts(mh.level_parts[1],GridapDistributed.local_views(uhold),
    #                  GridapDistributed.local_views(dΩ_old)) do part,uhold,dΩ_old
    #    if part==1
    #      println("PPP $(get_cell_dof_values(uhold))")
    #    end
    #  end
    # map_parts(mh.level_parts[1],GridapDistributed.local_views(uhnew),
    #                  GridapDistributed.local_views(dΩ_new)) do part,uhnew,dΩ_new
    #    if part==2
    #      println("YYY $(get_cell_dof_values(uhnew))")
    #    end
    #  end
  end

  # Give me how many processors you want per level
  # in an array with as many entries as levels
  num_parts_x_level = [4,2,1]
  if !MPI.Initialized()
    MPI.Init()
  end
  parts = get_part_ids(mpi,4)
  run(parts,num_parts_x_level)
  MPI.Finalize()
end
