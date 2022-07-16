module GMGLinearSolverTests
  using MPI
  using Gridap
  using GridapPETSc
  using PartitionedArrays
  using GridapDistributed
  using GridapP4est
  using P4est_wrapper
  using Test
  using LinearAlgebra
  using IterativeSolvers

  function generate_fe_spaces(mh)
    order=1
    reffe=ReferenceFE(lagrangian,Float64,order)
    test_spaces = Vector{FESpaceHierarchyLevel}(undef,num_levels(mh))
    trial_spaces = Vector{FESpaceHierarchyLevel}(undef,num_levels(mh))
    for i=1:num_levels(mh)
      model = get_level_model(mh,i)
      if (GridapP4est.i_am_in(model.parts))
         Vh = TestFESpace(get_level(mh,i),reffe,dirichlet_tags="boundary")
         Uh = TrialFESpace(u,Vh)
         test_spaces[i]  = Vh
         trial_spaces[i] = Uh
      end
    end
    test_spaces, trial_spaces
  end

  function assemble_mass_matrix(model,Vh,Uh,degree)
    Ω  = Triangulation(model.dmodel)
    dΩ = Measure(Ω,degree)
    a(u,v)=∫(v⋅u)dΩ
    assemble_matrix(a,Uh,Vh)
  end

  # function generate_mass_matrices(mh, fespaces)
  #   tests,trials=fespaces
  #   order=1
  #   mmatrices=Vector{PSparseMatrix}(undef,num_levels(mh))
  #   for i=1:num_levels(mh)
  #     model = get_level_model_before_redist(mh,i)
  #     if (GridapP4est.i_am_in(model.parts))
  #       Vh = get_level_fe_space_before_redist(tests[i])
  #       Uh = get_level_fe_space_before_redist(trials[i])
  #       A = assemble_mass_matrix(model,Vh,Uh,order)
  #     end
  #   end
  #   mmatrices
  # end

  function generate_stiffness_matrices(mh, fespaces)
    tests,trials=fespaces
    Gridap.Helpers.@check num_levels(mh)==length(tests)
    Gridap.Helpers.@check num_levels(mh)==length(trials)
    order=1
    matrices=Vector{PSparseMatrix}(undef,num_levels(mh))
    for i=1:num_levels(mh)
      model = get_level_model(mh,i)
      if (GridapP4est.i_am_in(model.parts))
        Ω  = Triangulation(model.dmodel)
        dΩ = Measure(Ω,2*(order+1))
        Vh = get_level_fe_space(tests[i])
        Uh = get_level_fe_space(trials[i])
        a(u,v)=∫(∇(v)⋅∇(u))dΩ
        A = assemble_matrix(a,Uh,Vh)
        matrices[i]=A
      end
    end
    matrices
  end

  mutable struct InterpolationMat
    Ωh_red
    dΩh_red
    Ωh
    ΩH
    dΩh
    UH
    VH
    Uh
    Vh
    Uh_red
    Vh_red
    mesh_hierarchy_level
    mesh_hierarchy_level_next
    uH_zero_dirichlet_values
    dof_values_H_fe_space_layout
    uh_zero_dirichlet_values
    dof_values_h_fe_space_layout
    dof_values_h_sys_layout_b
    dof_values_h_sys_layout_x
    mh
  end

  # TO-DO: InterpolationMatCache
  function LinearAlgebra.mul!(x::PVector,A::InterpolationMat,y::Union{PVector,Nothing})
    parts = A.mesh_hierarchy_level.model.parts
    Gridap.Helpers.@check GridapP4est.i_am_in(parts)

    model_next=get_level_model(A.mesh_hierarchy_level_next)
    coarse_parts  = model_next.parts
    i_am_in_coarse= GridapP4est.i_am_in(coarse_parts)

    if (i_am_in_coarse)
      # To-think: Do we need communication after copy!?
      copy!(A.dof_values_H_fe_space_layout,y)
      exchange!(A.dof_values_H_fe_space_layout)
      uH = FEFunction(A.UH,
                      A.dof_values_H_fe_space_layout,
                      A.uH_zero_dirichlet_values)
      # writevtk(A.ΩH,"uH",nsubcells=2,cellfields=["uH"=>uH])
    else
      uH = nothing
    end

    uH_h = change_domain_coarse_to_fine(A.Uh,uH,A.Ωh,A.mesh_hierarchy_level.ref_glue)
    #writevtk(A.Ωh,"uH_h",cellfields=["uH_h"=>uH_h])

    l(v) = ∫(v*uH_h)A.dΩh
    Gridap.FESpaces.assemble_vector!(l,A.dof_values_h_sys_layout_b,A.Vh)
    parts=A.mesh_hierarchy_level.model.parts

    if (A.mesh_hierarchy_level.model_red != nothing)
      fill!(A.dof_values_h_sys_layout_x,0.0)
      # println("aaaaaaaaaa aaaaaaaaaaaaaa")
      IterativeSolvers.cg!(A.dof_values_h_sys_layout_x,
                           A.mh,
                           A.dof_values_h_sys_layout_b;
                           verbose=i_am_main(parts),
                           reltol=1.0e-14)
      # println("bbbbbbbbbb bbbbbbbbbbbbbb")
      copy!(A.dof_values_h_fe_space_layout,
            A.dof_values_h_sys_layout_x)
      exchange!(A.dof_values_h_fe_space_layout)

      # map_parts(A.dof_values_h_fe_space_layout.values,
      #           A.dof_values_h_sys_layout_x.values) do a,b
      #   println(a)
      #   println(b)
      # end

      uhold=FEFunction(A.Uh,
                       A.dof_values_h_fe_space_layout,
                       A.uh_zero_dirichlet_values)

      # println("XXX $(sum(∫(uhold)A.dΩh))")

      # map_parts(A.mesh_hierarchy_level.model.parts,uhold.fields) do part,fs
      #     println("AAA1 $(fs.cell_dof_values)")
      #     println("AAA1 $(fs.free_values)")
      #     println("AAA1 $(fs.dirichlet_values)")
      #     println("AAA1 $(fs.fe_space)")
      # end

      # writevtk(A.Ωh,"Correction",cellfields=["correction"=>uhold])

      uhnew  = redistribute_fe_function(uhold,
                                        A.Uh,
                                        A.Uh_red,
                                        A.mesh_hierarchy_level.model_red.dmodel,
                                        A.mesh_hierarchy_level.red_glue)

      #writevtk(A.Ωh_red,"Correction red",cellfields=["correction red"=>uhnew])

      # map_parts(A.mesh_hierarchy_level.model_red.parts,uhnew.fields) do part,fs
      #   if part == 1
      #     println("AAA $(fs.cell_dof_values)")
      #     println("AAA $(fs.free_values)")
      #     println("AAA $(fs.dirichlet_values)")
      #     println("AAA $(fs.fe_space)")
      #   end
      # end

      # println("YYY $(sum(∫(uhnew)A.dΩh_red))")


      # map_parts(uhnew.metadata.free_values.values) do fv
      #   println("RRR $(fv)")
      # end



      # COMM after copy???
      copy!(x,uhnew.metadata.free_values)
      exchange!(x)
      # map_parts(x.values,uhnew.metadata.free_values.values) do x1, x2
      #   println("XXX $(x1)")
      #   println("YYY $(x2)")
      # end
    else
      IterativeSolvers.cg!(x,
                           A.mh,
                           A.dof_values_h_sys_layout_b;
                           verbose=i_am_main(parts),
                           reltol=1.0e-14)
    end
  end

  function setup_interpolation_mat(mh,fespaces,level)
    Gridap.Helpers.@check 1 <= level <= num_levels(mh)-1

    tests,trials=fespaces
    order=1

    # Fine level (all processes participate)
    model_h=get_level_model_before_redist(mh,level)
    Ωh=Triangulation(model_h.dmodel)
    Ωh_red=Triangulation(get_level_model(mh,level).dmodel)
    dΩh=Measure(Ωh,2*(order+1))
    dΩh_red=Measure(Ωh_red,2*(order+1))


    Vh=get_level_fe_space_before_redist(tests[level])
    Uh=get_level_fe_space_before_redist(trials[level])
    Mh=assemble_mass_matrix(model_h,Vh,Uh,2*order+1)


    dof_values_h_fe_space_layout = PVector(0.0,Vh.gids)
    dof_values_h_sys_layout_b = similar(dof_values_h_fe_space_layout,(axes(Mh)[1],))
    dof_values_h_sys_layout_x = similar(dof_values_h_fe_space_layout,(axes(Mh)[2],))
    uh_zero_dirichlet_values=
       map_parts(GridapDistributed.local_views(Uh.spaces)) do space
          zeros(num_dirichlet_dofs(space))
       end

    # Coarse level (only processes of the next level participate)
    model_H=get_level_model(mh,level+1)
    if (GridapP4est.i_am_in(model_H.parts))
      ΩH = Triangulation(model_H.dmodel)
      VH=get_level_fe_space(tests[level+1])
      UH=get_level_fe_space(trials[level+1])
      dof_values_H_fe_space_layout = PVector(0.0,VH.gids)
      uH_zero_dirichlet_values=
         map_parts(GridapDistributed.local_views(UH.spaces)) do space
             zeros(num_dirichlet_dofs(space))
         end
      map_parts(mh.levels[level].model.parts) do part_id
          println("INTERP $(part_id) $(ΩH)")
      end
    else
      ΩH = nothing
      UH = nothing
      VH = nothing
      dof_values_H_fe_space_layout = nothing
      uH_zero_dirichlet_values = nothing
      map_parts(mh.levels[level].model.parts) do part_id
        println("INTERP $(part_id) $(ΩH)")
      end
    end

    cache=InterpolationMat(Ωh_red,
                           dΩh_red,
                           Ωh,
                           ΩH,
                           dΩh,
                           UH,
                           VH,
                           Uh,
                           Vh,
                           get_level_fe_space(trials[level]),
                           get_level_fe_space(tests[level]),
                           get_level(mh,level),
                           get_level(mh,level+1),
                           uH_zero_dirichlet_values,
                           dof_values_H_fe_space_layout,
                           uh_zero_dirichlet_values,
                           dof_values_h_fe_space_layout,
                           dof_values_h_sys_layout_b,
                           dof_values_h_sys_layout_x,
                           Mh)
  end

  mutable struct RestrictionMat
    Ωh
    ΩH
    dΩh
    dΩH
    Uh
    Uh_red
    UH
    VH
    mesh_hierarchy_level
    mesh_hierarchy_level_next
    uh_zero_dirichlet_values_red
    dof_values_h_fe_space_layout_red
    uH_zero_dirichlet_values
    dof_values_H_sys_layout_b
    mH
  end

  function LinearAlgebra.mul!(x::Union{PVector,Nothing}, A::RestrictionMat, y::PVector)
    parts = A.mesh_hierarchy_level.model.parts
    Gridap.Helpers.@check GridapP4est.i_am_in(parts)

    copy!(A.dof_values_h_fe_space_layout_red,y)
    exchange!(A.dof_values_h_fe_space_layout_red)
    if (A.mesh_hierarchy_level.model_red != nothing)
      uhold = FEFunction(A.Uh_red,
                         A.dof_values_h_fe_space_layout_red,
                         A.uh_zero_dirichlet_values_red)

      # x is in a subcommunicator of y
      uh=redistribute_fe_function(uhold,
                                  A.Uh_red,
                                  A.Uh,
                                  A.mesh_hierarchy_level.model.dmodel,
                                  A.mesh_hierarchy_level.red_glue;
                                  reverse=true)

      # writevtk(A.Ωh,"Residual",cellfields=["residual"=>uh])
    else
      # x and y are in the same communicator
      uh = FEFunction(A.Uh,
                      A.dof_values_h_fe_space_layout_red,
                      A.uh_zero_dirichlet_values_red)
    end

    # uh is in h communicator, but with void parts for those tasks not in the next level
    # uh_H is void for those tasks not in the next level
    uh_H = change_domain_fine_to_coarse(uh,A.ΩH,A.mesh_hierarchy_level.ref_glue)

    parts = A.mesh_hierarchy_level_next.model.parts
    if (GridapP4est.i_am_in(parts))
      Gridap.Helpers.@check x != nothing
      l(v) = ∫(v*uh_H)A.dΩH
      Gridap.FESpaces.assemble_vector!(l,A.dof_values_H_sys_layout_b,A.VH)
      fill!(x,0.0)
      IterativeSolvers.cg!(x,
                           A.mH,
                           A.dof_values_H_sys_layout_b;
                           verbose=i_am_main(parts),
                           reltol=1.0e-14)
      uH = FEFunction(A.UH,
                      x,
                      A.uH_zero_dirichlet_values)
      # writevtk(A.ΩH,"Residual_projected",cellfields=["residual_projected"=>uH])
    end
  end

  function setup_restriction_mat(mh,fespaces,level)
    Gridap.Helpers.@check 1 <= level <= num_levels(mh)-1
    tests,trials=fespaces
    order=1

    model_h=get_level_model_before_redist(mh,level)
    Ωh  = Triangulation(model_h.dmodel)
    dΩh = Measure(Ωh,2*order+1)
    Uh  = get_level_fe_space_before_redist(trials[level])
    Vh  = get_level_fe_space_before_redist(tests[level])
    Uh_red= get_level_fe_space(trials[level])
    Vh_red= get_level_fe_space(tests[level])
    dof_values_h_fe_space_layout_red = PVector(0.0,Vh_red.gids)
    dof_values_h_fe_space_layout     = PVector(0.0,Vh.gids)
    uh_zero_dirichlet_values_red=
      map_parts(GridapDistributed.local_views(Uh_red.spaces)) do space
           zeros(num_dirichlet_dofs(space))
      end

    model_H=get_level_model(mh,level+1)
    if (GridapP4est.i_am_in(model_H.parts))
      ΩH  = Triangulation(model_H.dmodel)
      dΩH = Measure(ΩH,2*(order+1))
      VH=get_level_fe_space(tests[level+1])
      UH=get_level_fe_space(trials[level+1])
      MH=assemble_mass_matrix(model_H,VH,UH,2*order+1)
      # Generate arrays of DoF Values for coarse system
      dof_values_H_fe_space_layout = PVector(0.0,VH.gids)
      dof_values_H_sys_layout_b = similar(dof_values_H_fe_space_layout,(axes(MH)[1],))
      uH_zero_dirichlet_values=
      map_parts(GridapDistributed.local_views(UH.spaces)) do space
           zeros(num_dirichlet_dofs(space))
      end
      map_parts(mh.levels[level].model.parts) do part_id
        println("RESTRICT $(part_id) $(ΩH)")
      end
    else
      ΩH  = nothing
      dΩH = nothing
      VH  = nothing
      UH  = nothing
      MH  = nothing
      dof_values_H_fe_space_layout = nothing
      dof_values_H_sys_layout_b    = nothing
      uH_zero_dirichlet_values = nothing
      map_parts(mh.levels[level].model.parts) do part_id
        println("RESTRICT $(part_id) $(ΩH)")
      end
    end

    cache=RestrictionMat(Ωh,
                         ΩH,
                         dΩh,
                         dΩH,
                         Uh,
                         Uh_red,
                         UH,
                         VH,
                         get_level(mh,level),
                         get_level(mh,level+1),
                         uh_zero_dirichlet_values_red,
                         dof_values_h_fe_space_layout_red,
                         uH_zero_dirichlet_values,
                         dof_values_H_sys_layout_b,
                         MH)
  end

  function setup_interpolations_and_restrictions(mh,fespaces)
    nlevs = num_levels(mh)
    interpolations=Vector{InterpolationMat}(undef,nlevs-1)
    restrictions=Vector{RestrictionMat}(undef,nlevs-1)
    for l=1:nlevs-1
      model = get_level_model(mh,l)
      if (GridapP4est.i_am_in(model.parts))
        interpolations[l]=setup_interpolation_mat(mh,fespaces,l)
        restrictions[l]=setup_restriction_mat(mh,fespaces,l)
      end
    end
    interpolations, restrictions
  end

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

    num_parts_x_level = [2,1]
    cmodel=CartesianDiscreteModel(domain,(4,4))
    mh=ModelHierarchy(parts,cmodel,num_parts_x_level)
    fespaces=generate_fe_spaces(mh)
    smatrices=generate_stiffness_matrices(mh,fespaces)
    interp,restrict=setup_interpolations_and_restrictions(mh,fespaces)

    model=get_level_model(mh,1)
    Ω = Triangulation(model.dmodel)
    order=1
    dΩ = Measure(Ω,2*(order+1))
    tests,trials=fespaces
    Vh = tests[1]
    Uh = trials[1]
    a(u,v)=∫(∇(v)⋅∇(u))dΩ
    l(v)=∫(v*f)dΩ
    Vh = get_level_fe_space(tests[1])
    Uh = get_level_fe_space(trials[1])

    op=AffineFEOperator(a,l,Uh,Vh)
    A=op.op.matrix
    b=op.op.vector
    x=PVector(0.0,A.cols)

    GMG!(x,
         b,
         mh,
         fespaces,
         smatrices,
         interp,
         restrict;
         rtol=1.0e-06,
         maxiter=200,
         smooth_iter=5)

    uh=FEFunction(Uh,x)
    # Error norms and print solution
    e = u-uh
    e_l2 = sum(∫(e*e)dΩ)
    tol = 1.0e-9
    @test e_l2 < tol
    map_parts(parts) do part
    if (part==1)
      println("$(e_l2) < $(tol)\n")
    end
    end
    model_hierarchy_free!(mh)
  end

  if !MPI.Initialized()
    MPI.Init()
  end
  parts = get_part_ids(mpi,2)
  run(parts,(1,1))
  MPI.Finalize()
end
