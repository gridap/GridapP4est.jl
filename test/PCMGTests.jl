module PCMGTests
  using MPI
  using Gridap
  using GridapPETSc
  using PartitionedArrays
  using GridapDistributed
  using GridapP4est
  using P4est_wrapper
  using Test

  function generate_fe_spaces(mh)
    order=1
    reffe=ReferenceFE(lagrangian,Float64,order)
    test_spaces = Vector{GridapP4est.FESpaceHierarchyLevel}(undef,num_levels(mh))
    trial_spaces = Vector{GridapP4est.FESpaceHierarchyLevel}(undef,num_levels(mh))
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

  function generate_mass_matrices(model_hierarchy, fespaces)
    tests,trials=fespaces
    order=1
    # Here we are assuming that the first level of the
    # hierarchy is distributed among all processes
    model=get_level_model_before_redist(model_hierarchy,1)
    Ω  = Triangulation(model.dmodel)
    dΩ = Measure(Ω,2*(order+1))
    Vh = get_level_fe_space_before_redist(tests[1])
    Uh = get_level_fe_space_before_redist(trials[1])
    a(u,v)=∫(v*u)dΩ
    A = assemble_matrix(a,Uh,Vh)
    mmatrices=Vector{typeof(A)}(undef,num_levels(model_hierarchy))
    mmatrices[1]=A
    for i=2:num_levels(model_hierarchy)
      model = get_level_model_before_redist(model_hierarchy,i)
      if (GridapP4est.i_am_in(model.parts))
        Ω  = Triangulation(model.dmodel)
        dΩ = Measure(Ω,2*(order+1))
        Vh = get_level_fe_space_before_redist(tests[i])
        Uh = get_level_fe_space_before_redist(trials[i])
        a(u,v)=∫(v*u)dΩ
        A = assemble_matrix(a,Uh,Vh)
        mmatrices[i]=A
      end
    end
    mmatrices
  end

  function generate_stiffness_matrices(model_hierarchy, fespaces)
    tests,trials=fespaces
    Gridap.Helpers.@check num_levels(model_hierarchy)==length(tests)
    Gridap.Helpers.@check num_levels(model_hierarchy)==length(trials)
    order=1
    matrices=Vector{PETScMatrix}(undef,num_levels(model_hierarchy))
    for i=1:num_levels(model_hierarchy)
      model = get_level_model(model_hierarchy,i)
      if (GridapP4est.i_am_in(model.parts))
        Ω  = Triangulation(model.dmodel)
        dΩ = Measure(Ω,2*(order+1))
        Vh = get_level_fe_space(tests[i])
        Uh = get_level_fe_space(trials[i])
        a(u,v)=∫(∇(v)⋅∇(u))dΩ
        A = assemble_matrix(a,Uh,Vh)
        A = convert(PETScMatrix,A)
        matrices[i]=A
      end
    end
    matrices
  end

  mutable struct InterpolationMatCache
    Ωh
    ΩH
    dΩh
    UH
    VH
    Uh
    Vh
    hierarchy_level
    uH_zero_dirichlet_values
    dof_values_H_fe_space_layout
    dof_values_H_sys_layout
    dof_values_h_sys_layout_b
    dof_values_h_sys_layout_x
    ns
  end

  function matmult_interpolation(A::Ptr{Cvoid}, x::Ptr{Cvoid}, y::Ptr{Cvoid})
    ctx = Ref{Ptr{Cvoid}}()
    @check_error_code GridapP4est.MatShellGetContext(A,ctx)
    cache  = unsafe_pointer_to_objref(ctx[])

    # @check_error_code GridapPETSc.PETSC.VecView(Vec(x),C_NULL)

    GridapPETSc._copy!(cache.dof_values_H_sys_layout, Vec(x))

    # map_parts(cache.dof_values_H_sys_layout.values) do u
    #   println(u)
    # end

    GridapPETSc.copy!(cache.dof_values_H_fe_space_layout,cache.dof_values_H_sys_layout)

    # map_parts(cache.dof_values_H_sys_layout.values) do u
    #   println(u)
    # end

    uH = FEFunction(cache.UH,
                    cache.dof_values_H_fe_space_layout,
                    cache.uH_zero_dirichlet_values)

    uH_h = change_domain_coarse_to_fine(uH,cache.Ωh,cache.hierarchy_level.ref_glue)

    #writevtk(cache.ΩH,"uH",nsubcells=2,cellfields=["uH"=>uH])
    #writevtk(cache.Ωh,"uH_h",cellfields=["uH_h"=>uH_h])

    l(v) = ∫(v*uH_h)cache.dΩh
    Gridap.FESpaces.assemble_vector!(l,cache.dof_values_h_sys_layout_b,cache.Vh)

    # map_parts(cache.dof_values_h_sys_layout_b.values) do u
    #   println(u)
    # end

    fill!(cache.dof_values_h_sys_layout_x,0.0)
    solve!(cache.dof_values_h_sys_layout_x,cache.ns,cache.dof_values_h_sys_layout_b)

    # map_parts(cache.ns.A.values,cache.dof_values_h_sys_layout_b.values) do A,b
    #  println(A)
    #  println(b)
    #  println(A\b)
    # end

    # uh = FEFunction(cache.Uh,cache.dof_values_h_sys_layout_x)
    # writevtk(cache.Ωh,"uh",cellfields=["uh"=>uh])

    GridapPETSc._copy!(Vec(y),cache.dof_values_h_sys_layout_x)
    # println("MatMult Interp")
    # @check_error_code GridapPETSc.PETSC.VecView(Vec(y),C_NULL)

    GridapPETSc.PETSC.PetscErrorCode(0)
  end

  function set_ksp_mm(ksp)
    @check_error_code GridapPETSc.PETSC.KSPSetOptionsPrefix(ksp[],"intergrid_mm_")
    @check_error_code GridapPETSc.PETSC.KSPSetFromOptions(ksp[])
    #@check_error_code GridapPETSc.PETSC.KSPView(ksp[],C_NULL)
  end

  function setup_interpolation_mat(mh,fespaces,mmatrices,level)
    Gridap.Helpers.@check 1 <= level <= num_levels(mh)-1

    tests,trials=fespaces
    order=1

    # Fine level (all processes participate)
    model_h=get_level_model_before_redist(mh,level)
    Ωh=Triangulation(model_h.dmodel)
    dΩh=Measure(Ωh,2*(order+1))
    Vh=get_level_fe_space_before_redist(tests[level])
    Uh=get_level_fe_space_before_redist(trials[level])
    Mh=mmatrices[level]
    dof_values_h_fe_space_layout = PVector(0.0,Vh.gids)
    dof_values_h_sys_layout_b = similar(dof_values_h_fe_space_layout,(axes(Mh)[2],))
    dof_values_h_sys_layout_x = similar(dof_values_h_fe_space_layout,(axes(Mh)[2],))

    # Coarse level (only processes of the next level participate)
    model_H=get_level_model(mh,level+1)
    if (GridapP4est.i_am_in(model_H.parts))
      ΩH = Triangulation(model_H.dmodel)
      VH=get_level_fe_space(tests[level+1])
      UH=get_level_fe_space(trials[level+1])
      MH=mmatrices[level+1]
      dof_values_H_fe_space_layout = PVector(0.0,VH.gids)
      dof_values_H_sys_layout = similar(dof_values_H_fe_space_layout,(axes(MH)[2],))
      uH_zero_dirichlet_values=
         map_parts(GridapDistributed.local_views(UH.spaces)) do space
             zeros(num_dirichlet_dofs(space))
         end
    else
      ΩH = nothing
      UH = nothing
      VH = nothing
      dof_values_H_fe_space_layout = nothing
      dof_values_H_sys_layout = nothing
      uH_zero_dirichlet_values = nothing
    end

    # To-think: perhaps we may pre-compute this
    solver=PETScLinearSolver(set_ksp_mm)
    ss=symbolic_setup(solver,Mh)
    ns=numerical_setup(ss,Mh)

    cache=InterpolationMatCache(Ωh,
                                ΩH,
                                dΩh,
                                UH,
                                VH,
                                Uh,
                                Vh,
                                get_level(mh,level),
                                uH_zero_dirichlet_values,
                                dof_values_H_fe_space_layout,
                                dof_values_H_sys_layout,
                                dof_values_h_sys_layout_b,
                                dof_values_h_sys_layout_x,
                                ns)

    ctx = pointer_from_objref(cache)
    fptr = @cfunction(matmult_interpolation,
                      GridapPETSc.PETSC.PetscErrorCode,
                      (Ptr{Cvoid},Ptr{Cvoid},Ptr{Cvoid}))

    comm = model_h.parts.comm
    petscmat=PETScMatrix(comm)

    m=0
    map_parts(Vh.gids.partition) do gids
      m = num_oids(gids)
    end

    n=0
    if (GridapP4est.i_am_in(model_H.parts))
      map_parts(VH.gids.partition) do gids
        n = num_oids(gids)
      end
    end

    println("INTERPOLATION m=$(m) n=$(n)")

    M = GridapPETSc.PETSC.PETSC_DETERMINE
    N = GridapPETSc.PETSC.PETSC_DETERMINE
    @check_error_code GridapP4est.MatCreateShell(comm,m,n,M,N,ctx,petscmat.mat)
    @check_error_code GridapP4est.MatShellSetOperation(petscmat.mat[],GridapP4est.MATOP_MULT,fptr)
    petscmat.ownership = cache
    petscmat
  end

  mutable struct RestrictionMatCache
    Ωh
    ΩH
    dΩH
    Uh
    UH
    VH
    hierarchy_level
    uh_zero_dirichlet_values
    dof_values_h_fe_space_layout
    dof_values_H_sys_layout_b
    dof_values_H_sys_layout_x
    dof_values_h_sys_layout
    ns
  end

  function matmult_restriction(A::Ptr{Cvoid}, x::Ptr{Cvoid}, y::Ptr{Cvoid})
    ctx = Ref{Ptr{Cvoid}}()
    @check_error_code GridapP4est.MatShellGetContext(A,ctx)
    cache  = unsafe_pointer_to_objref(ctx[])
    # @check_error_code GridapPETSc.PETSC.VecView(Vec(x),C_NULL)
    GridapPETSc._copy!(cache.dof_values_h_sys_layout, Vec(x))

    # map_parts(cache.dof_values_h_sys_layout.values) do u
    #   println(u)
    # end

    GridapPETSc.copy!(cache.dof_values_h_fe_space_layout,cache.dof_values_h_sys_layout)

    # map_parts(cache.dof_values_h_fe_space_layout.values) do u
    #   println(u)
    # end

    uh = FEFunction(cache.Uh,
                    cache.dof_values_h_fe_space_layout,
                    cache.uh_zero_dirichlet_values)

    uh_H = change_domain_fine_to_coarse(uh,cache.ΩH,cache.hierarchy_level.ref_glue)

    #writevtk(cache.Ωh,"uh",cellfields=["uh"=>uh])
    #writevtk(cache.ΩH,"uh_H",nsubcells=2,cellfields=["uh_H"=>uh_H])

    l(v) = ∫(v*uh_H)cache.dΩH
    Gridap.FESpaces.assemble_vector!(l,cache.dof_values_H_sys_layout_b,cache.VH)

    # map_parts(cache.dof_values_H_sys_layout_b.values) do u
    #   println(u)
    # end

    fill!(cache.dof_values_H_sys_layout_x,0.0)
    solve!(cache.dof_values_H_sys_layout_x,cache.ns,cache.dof_values_H_sys_layout_b)

    # map_parts(cache.ns.A.values,cache.dof_values_H_sys_layout_b.values) do A,b
    #   println(A)
    #   println(b)
    #   println(A\b)
    #  end

    # map_parts(cache.dof_values_H_sys_layout_x.values) do u
    #   println(u)
    # end

    # uH = FEFunction(cache.UH,cache.dof_values_H_sys_layout_x)
    #writevtk(cache.ΩH,"uH",cellfields=["uH"=>uH])

    GridapPETSc._copy!(Vec(y),cache.dof_values_H_sys_layout_x)
    # @check_error_code GridapPETSc.PETSC.VecView(Vec(y),C_NULL)

    GridapPETSc.PETSC.PetscErrorCode(0)
  end

  function setup_restriction_mat(mh,fespaces,mmatrices,level)
    Gridap.Helpers.@check 1 <= level <= num_levels(mh)-1
    tests,trials=fespaces
    order=1

    model_h=get_level_model_before_redist(mh,level)
    Ωh=Triangulation(model_h.dmodel)
    Mh=mmatrices[level]
    Vh=get_level_fe_space_before_redist(tests[level])
    Uh=get_level_fe_space_before_redist(trials[level])
    dof_values_h_fe_space_layout = PVector(0.0,Vh.gids)
    dof_values_h_sys_layout = similar(dof_values_h_fe_space_layout,(axes(Mh)[2],))
    uh_zero_dirichlet_values=
      map_parts(GridapDistributed.local_views(Uh.spaces)) do space
           zeros(num_dirichlet_dofs(space))
      end

    model_H=get_level_model(mh,level+1)
    if (GridapP4est.i_am_in(model_H.parts))
      ΩH  = Triangulation(model_H.dmodel)
      dΩH = Measure(ΩH,2*(order+1))
      VH=get_level_fe_space(tests[level+1])
      UH=get_level_fe_space(trials[level+1])
      MH=mmatrices[level+1]
      # Generate arrays of DoF Values for fine system
      dof_values_H_fe_space_layout = PVector(0.0,VH.gids)
      dof_values_H_sys_layout_b = similar(dof_values_H_fe_space_layout,(axes(MH)[2],))
      dof_values_H_sys_layout_x = similar(dof_values_H_fe_space_layout,(axes(MH)[2],))

      # To-think: perhaps we may pre-compute this
      solver=PETScLinearSolver(set_ksp_mm)
      ss=symbolic_setup(solver,MH)
      ns=numerical_setup(ss,MH)
    else
      ΩH  = nothing
      dΩH = nothing
      VH  = nothing
      UH  = nothing
      MH  = nothing
      dof_values_H_fe_space_layout = nothing
      dof_values_H_sys_layout_b    = nothing
      dof_values_H_sys_layout_x    = nothing
      ns  = nothing
    end

    cache=RestrictionMatCache(Ωh,
                              ΩH,
                              dΩH,
                              Uh,
                              UH,
                              VH,
                              get_level(mh,level),
                              uh_zero_dirichlet_values,
                              dof_values_h_fe_space_layout,
                              dof_values_H_sys_layout_b,
                              dof_values_H_sys_layout_x,
                              dof_values_h_sys_layout,
                              ns)

    ctx = pointer_from_objref(cache)
    fptr = @cfunction(matmult_restriction,
                      GridapPETSc.PETSC.PetscErrorCode,
                      (Ptr{Cvoid},Ptr{Cvoid},Ptr{Cvoid}))

    n=0
    map_parts(Vh.gids.partition) do gids
      n = num_oids(gids)
    end

    m=0
    if (GridapP4est.i_am_in(model_H.parts))
      map_parts(VH.gids.partition) do gids
        m = num_oids(gids)
      end
    end

    println("RESTRICTION m=$(m) n=$(n)")

    comm = model_h.parts.comm
    petscmat=PETScMatrix(comm)
    M = GridapPETSc.PETSC.PETSC_DETERMINE
    N = GridapPETSc.PETSC.PETSC_DETERMINE
    @check_error_code GridapP4est.MatCreateShell(comm,m,n,M,N,ctx,petscmat.mat)
    @check_error_code GridapP4est.MatShellSetOperation(petscmat.mat[],GridapP4est.MATOP_MULT,fptr)
    petscmat.ownership = cache
    petscmat
  end

  function setup_interpolations_and_restrictions(mh,fespaces,mmatrices)
    nlevs = num_levels(mh)
    interpolations=Vector{PETScMatrix}(undef,nlevs-1)
    restrictions=Vector{PETScMatrix}(undef,nlevs-1)
    for l=1:nlevs-1
      model = get_level_model(mh,l)
      if (GridapP4est.i_am_in(model.parts))
        interpolations[l]=setup_interpolation_mat(mh,fespaces,mmatrices,l)
        restrictions[l]=setup_restriction_mat(mh,fespaces,mmatrices,l)
      end
    end
    interpolations, restrictions
  end

  function setup_KSP_PCMG!(ksp,
                           mh,
                           fespaces,
                           mmatrices,
                           smatrices,
                           interpolations,
                           restrictions)
    nlevs=num_levels(mh)
    println("NLs $(nlevs)")
    pc=Ref{PC}()
    ksp_smoother=Ref{KSP}()
    @check_error_code GridapPETSc.PETSC.KSPSetFromOptions(ksp[])
    @check_error_code GridapPETSc.PETSC.KSPGetPC(ksp[],pc)
    @check_error_code GridapPETSc.PETSC.PCSetType(pc[],GridapPETSc.PETSC.PCMG)


    comms = [ mh.level_parts[i].comm.val for i=nlevs:-1:1 ]
    @check_error_code GridapP4est.PCMGSetLevels(pc[],PetscInt(nlevs),comms)

    for l=1:nlevs-1
      petsc_level=nlevs-l
      model=get_level_model(mh,l)
      if (GridapP4est.i_am_in(model.parts))
        @check_error_code GridapP4est.PCMGSetInterpolation(pc[],
                                            PetscInt(petsc_level),
                                            interpolations[l].mat[])
        @check_error_code GridapP4est.PCMGSetRestriction(pc[],
                                            PetscInt(petsc_level),
                                            restrictions[l].mat[])
        if l==1
          @check_error_code GridapPETSc.PETSC.KSPSetOperators(ksp[],
                              smatrices[l].mat[],
                              smatrices[l].mat[])
        end

        @check_error_code GridapP4est.PCMGGetSmoother(pc[],
                                                      PetscInt(petsc_level),
                                                      ksp_smoother)

        @check_error_code GridapPETSc.PETSC.KSPSetOperators(ksp_smoother[],
                                                            smatrices[l].mat[],
                                                            smatrices[l].mat[])
      end
    end
    model=get_level_model(mh,nlevs)
    if (GridapP4est.i_am_in(model.parts))
      @check_error_code GridapP4est.PCMGGetCoarseSolve(pc[],ksp_smoother);
      @check_error_code GridapPETSc.PETSC.KSPSetOperators(ksp_smoother[],
                                                          smatrices[nlevs].mat[],
                                                          smatrices[nlevs].mat[])
    else
      @check_error_code GridapPETSc.PETSC.KSPSetOperators(ksp_smoother[],
                                                    smatrices[1].mat[],
                                                    smatrices[1].mat[])
    end
    @check_error_code GridapPETSc.PETSC.KSPSetUp(ksp[])
    @check_error_code GridapPETSc.PETSC.KSPView(ksp[],C_NULL)
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

    options = "-ksp_type cg -ksp_monitor -ksp_converged_reason -ksp_rtol 1.0e-8 -pc_type mg -pc_mg_type multiplicative -pc_mg_multiplicative_cycles 1  -mg_coarse_pc_type cholesky -mg_coarse_ksp_type preonly -mg_levels_pc_type jacobi -mg_levels_ksp_type richardson -intergrid_mm_ksp_type preonly -intergrid_mm_pc_type cholesky"

    num_parts_x_level = [2,1]
    GridapPETSc.with(args=split(options)) do
      cmodel=CartesianDiscreteModel(domain,(2,2))
      mh=ModelHierarchy(parts,cmodel,num_parts_x_level)
      fespaces=generate_fe_spaces(mh)
      mmatrices=generate_mass_matrices(mh,fespaces)
      smatrices=generate_stiffness_matrices(mh,fespaces)
      interpolations,restrictions=setup_interpolations_and_restrictions(mh,
                                                                        fespaces,
                                                                        mmatrices)
      # Setup solver via low level PETSC API calls
      function mykspsetup(ksp)
        setup_KSP_PCMG!(ksp,
                        mh,
                        fespaces,
                        mmatrices,
                        smatrices,
                        interpolations,
                        restrictions)
      end
      ls = PETScLinearSolver(mykspsetup)
      fels = LinearFESolver(ls)

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
      uh = solve(fels,op)

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
  end
  if !MPI.Initialized()
    MPI.Init()
  end
  parts = get_part_ids(mpi,2)
  run(parts,(1,1))
  MPI.Finalize()
end
