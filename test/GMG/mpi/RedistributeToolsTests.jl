module RedistributeToolsTests
  using MPI
  using PartitionedArrays
  using Gridap
  using GridapDistributed
  using GridapP4est
  using P4est_wrapper
  using Test

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

    uhold = interpolate(u,UNEW)
    uhnew  = redistribute_fe_function(uhold,
                                      UNEW,
                                      UOLD,
                                      mh.levels[1].model.dmodel,
                                      mh.levels[1].red_glue;
                                      reverse=true)

    o=sum(∫(uhnew)dΩ_old)
    n=sum(∫(uhold)dΩ_new)
    @test o ≈ n
    model_hierarchy_free!(mh)

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
  num_parts_x_level = [4,2,2]
  ranks=num_parts_x_level[1]
  prun(run,mpi,ranks,num_parts_x_level)
  MPI.Finalize()
end
