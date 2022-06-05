module PCMGTests
  using MPI
  using Gridap
  using GridapPETSc
  using PartitionedArrays
  using GridapDistributed
  using GridapP4est
  using P4est_wrapper
  using Test

  # TO-DO: generalize (order, reffe, etc.)
  function generate_meshes(parts,cmodel,num_levels)
    model=OctreeDistributedDiscreteModel(parts,cmodel)
    meshes=Vector{typeof(model)}(undef,num_levels)
    meshes[num_levels]=model
    fmodel,glue=refine(model)
    meshes[num_levels-1]=fmodel
    glues = Vector{typeof(glue)}(undef,num_levels-1)
    glues[num_levels-1] = glue
    for i=num_levels-2:-1:1
      fmodel,glue=refine(fmodel)
      meshes[i]=fmodel
      glues[i]=glue
    end
    meshes, glues
  end

  # TO-DO: generalize (order, reffe, etc.)
  function generate_fe_spaces(meshes)
    order=1
    reffe=ReferenceFE(lagrangian,Float64,order)
    Vh = TestFESpace(meshes[1].dmodel,reffe,dirichlet_tags="boundary")
    Uh = TrialFESpace(u,Vh)
    test_spaces     = Vector{typeof(Vh)}(undef,length(meshes))
    test_spaces[1]  = Vh
    trial_spaces    = Vector{typeof(Uh)}(undef,length(meshes))
    trial_spaces[1] = Uh
    for i=2:length(meshes)
      Vh = TestFESpace(meshes[i].dmodel,reffe,dirichlet_tags="boundary")
      Uh = TrialFESpace(u,Vh)
      test_spaces[i]  = Vh
      trial_spaces[i] = Uh
    end
    test_spaces, trial_spaces
  end

  function generate_mass_matrices(meshes, fespaces)
    tests,trials=fespaces
    order=1
    Ω  = Triangulation(meshes[1].dmodel)
    dΩ = Measure(Ω,2*(order+1))
    Vh = tests[1]
    Uh = trials[1]
    a(u,v)=∫(v*u)dΩ
    A = assemble_matrix(a,Uh,Vh)
    mmatrices=Vector{typeof(A)}(undef,length(meshes))
    mmatrices[1]=A
    for i=2:length(meshes)
      Ω  = Triangulation(meshes[i].dmodel)
      dΩ = Measure(Ω,2*(order+1))
      Vh = tests[i]
      Uh = trials[i]
      a(u,v)=∫(v*u)dΩ
      A = assemble_matrix(a,Uh,Vh)
      mmatrices[i]=A
    end
    mmatrices
  end

  function generate_stiffness_matrices(meshes, fespaces)
    tests,trials=fespaces
    Gridap.Helpers.@check length(meshes)==length(tests)
    Gridap.Helpers.@check length(meshes)==length(trials)
    order=1
    matrices=Vector{PETScMatrix}(undef,length(meshes))
    for i=1:length(meshes)
      Ω  = Triangulation(meshes[i].dmodel)
      dΩ = Measure(Ω,2*(order+1))
      Vh = tests[i]
      Uh = trials[i]
      a(u,v)=∫(∇(v)⋅∇(u))dΩ
      A = assemble_matrix(a,Uh,Vh)
      A = convert(PETScMatrix,A)
      matrices[i]=A
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
    glue
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

    uH_h = change_domain_coarse_to_fine(uH,cache.Ωh,cache.glue)

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

  function setup_interpolation_mat(parts,meshes,glues,fespaces,mmatrices,level)
    Gridap.Helpers.@check 1 <= level <= length(meshes)-1

    tests,trials=fespaces
    order=1
    Ωh = Triangulation(meshes[level].dmodel)
    ΩH = Triangulation(meshes[level+1].dmodel)
    dΩh=Measure(Ωh,2*(order+1))

    # Extract Vh/VH and mass matrices
    Vh=tests[level]
    Uh=trials[level]
    VH=tests[level+1]
    UH=trials[level+1]
    Mh=mmatrices[level]
    MH=mmatrices[level+1]

    # Generate arrays of DoF Values for coarse system
    dof_values_H_fe_space_layout = PVector(0.0,VH.gids)
    dof_values_h_fe_space_layout = PVector(0.0,Vh.gids)
    dof_values_H_sys_layout = similar(dof_values_H_fe_space_layout,(axes(MH)[2],))
    dof_values_h_sys_layout_b = similar(dof_values_h_fe_space_layout,(axes(Mh)[2],))
    dof_values_h_sys_layout_x = similar(dof_values_h_fe_space_layout,(axes(Mh)[2],))

    uH_zero_dirichlet_values=
      map_parts(GridapDistributed.local_views(UH.spaces)) do space
           zeros(num_dirichlet_dofs(space))
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
                                glues[level],
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
    m = PetscInt(num_free_dofs(Vh))
    n = PetscInt(num_free_dofs(VH))
    M = m
    N = n

    petscmat=PETScMatrix(parts.comm)
    @check_error_code GridapP4est.MatCreateShell(parts.comm,m,n,M,N,ctx,petscmat.mat)
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
    glue
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

    uh_H = change_domain_fine_to_coarse(uh,cache.ΩH,cache.glue)

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

  function setup_restriction_mat(parts,meshes,glues,fespaces,mmatrices,level)
    Gridap.Helpers.@check 1 <= level <= length(meshes)-1
    tests,trials=fespaces
    order=1
    Ωh  = Triangulation(meshes[level].dmodel)
    ΩH  = Triangulation(meshes[level+1].dmodel)
    dΩH = Measure(ΩH,2*(order+1))

    # Extract Vh/VH and mass matrices
    Vh=tests[level]
    VH=tests[level+1]
    Uh=trials[level]
    UH=trials[level+1]
    Mh=mmatrices[level]
    MH=mmatrices[level+1]

    # Generate arrays of DoF Values for fine system
    dof_values_h_fe_space_layout = PVector(0.0,Vh.gids)
    dof_values_H_fe_space_layout = PVector(0.0,VH.gids)
    dof_values_h_sys_layout = similar(dof_values_h_fe_space_layout,(axes(Mh)[2],))
    dof_values_H_sys_layout_b = similar(dof_values_H_fe_space_layout,(axes(MH)[2],))
    dof_values_H_sys_layout_x = similar(dof_values_H_fe_space_layout,(axes(MH)[2],))


    uh_zero_dirichlet_values=
      map_parts(GridapDistributed.local_views(Uh.spaces)) do space
           zeros(num_dirichlet_dofs(space))
      end

    # To-think: perhaps we may pre-compute this
    solver=PETScLinearSolver(set_ksp_mm)
    ss=symbolic_setup(solver,MH)
    ns=numerical_setup(ss,MH)

    cache=RestrictionMatCache(Ωh,
                              ΩH,
                              dΩH,
                              Uh,
                              UH,
                              VH,
                              glues[level],
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
    m = PetscInt(num_free_dofs(VH))
    n = PetscInt(num_free_dofs(Vh))
    M = m
    N = n
    petscmat=PETScMatrix(parts.comm)
    @check_error_code GridapP4est.MatCreateShell(parts.comm,m,n,M,N,ctx,petscmat.mat)
    @check_error_code GridapP4est.MatShellSetOperation(petscmat.mat[],GridapP4est.MATOP_MULT,fptr)
    petscmat.ownership = cache
    petscmat
  end

  function setup_interpolations_and_restrictions(parts,meshes,glues,fespaces,mmatrices)
    comm = parts.comm
    num_levels = length(meshes)
    interpolations=Vector{PETScMatrix}(undef,num_levels-1)
    restrictions=Vector{PETScMatrix}(undef,num_levels-1)
    for l=1:num_levels-1
      interpolations[l]=setup_interpolation_mat(parts,meshes,glues,fespaces,mmatrices,l)
      restrictions[l]=setup_restriction_mat(parts,meshes,glues,fespaces,mmatrices,l)
    end
    interpolations, restrictions
  end

  function setup_KSP_PCMG!(ksp,
                           parts,
                           meshes,
                           glues,
                           fespaces,
                           mmatrices,
                           smatrices,
                           interpolations,
                           restrictions)
    num_levels=length(meshes)
    comm = parts.comm
    pc=Ref{PC}()
    ksp_smoother=Ref{KSP}()
    @check_error_code GridapPETSc.PETSC.KSPSetFromOptions(ksp[])
    @check_error_code GridapPETSc.PETSC.KSPGetPC(ksp[],pc)
    @check_error_code GridapPETSc.PETSC.PCSetType(pc[],GridapPETSc.PETSC.PCMG)
    @check_error_code GridapP4est.PCMGSetLevels(pc[],PetscInt(num_levels),C_NULL)

    for l=1:num_levels-1
      petsc_level=num_levels-l
      @check_error_code GridapP4est.PCMGSetInterpolation(pc[],
                                           PetscInt(petsc_level),
                                           interpolations[l].mat[])
      @check_error_code GridapP4est.PCMGSetRestriction(pc[],
                                          PetscInt(petsc_level),
                                          restrictions[l].mat[])
      @check_error_code GridapP4est.PCMGGetSmoother(pc[],
                                      PetscInt(petsc_level),
                                      ksp_smoother)

      if l==1
        @check_error_code GridapPETSc.PETSC.KSPSetOperators(ksp[],
                             smatrices[l].mat[],
                             smatrices[l].mat[])
      end
      @check_error_code GridapPETSc.PETSC.KSPSetOperators(ksp_smoother[],
                                                            smatrices[l].mat[],
                                                            smatrices[l].mat[])
    end
    @check_error_code GridapP4est.PCMGGetCoarseSolve(pc[],ksp_smoother);
    @check_error_code GridapPETSc.PETSC.KSPSetOperators(ksp_smoother[],
                                                        smatrices[num_levels].mat[],
                                                        smatrices[num_levels].mat[])
    @check_error_code GridapPETSc.PETSC.KSPSetUp(ksp[])
    #@check_error_code GridapPETSc.PETSC.KSPView(ksp[],C_NULL)
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

    num_levels=2
    GridapPETSc.with(args=split(options)) do
       coarse_discrete_model=CartesianDiscreteModel(domain,(64,64))
       meshes,glues=generate_meshes(parts,coarse_discrete_model,num_levels)
       fespaces=generate_fe_spaces(meshes)
       mmatrices=generate_mass_matrices(meshes,fespaces)
       smatrices=generate_stiffness_matrices(meshes,fespaces)
       interpolations,restrictions=setup_interpolations_and_restrictions(parts,
                                                                         meshes,
                                                                         glues,
                                                                         fespaces,
                                                                         mmatrices)
       # Setup solver via low level PETSC API calls
       function mykspsetup(ksp)
          setup_KSP_PCMG!(ksp,
                          parts,
                          meshes,
                          glues,
                          fespaces,
                          mmatrices,
                          smatrices,
                          interpolations,
                          restrictions)
       end


       ls = PETScLinearSolver(mykspsetup)
       fels = LinearFESolver(ls)

       Ω    = Triangulation(meshes[1].dmodel)
       order=1
       dΩ = Measure(Ω,2*(order+1))
       tests,trials=fespaces
       Vh = tests[1]
       Uh = trials[1]
       a(u,v)=∫(∇(v)⋅∇(u))dΩ
       l(v)=∫(v*f)dΩ
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
