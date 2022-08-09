module GMGLinearSolverHDivRTTests
  using MPI
  using Gridap
  using Gridap.FESpaces
  using Gridap.ReferenceFEs
  using PartitionedArrays
  using GridapDistributed
  using GridapP4est
  using P4est_wrapper
  using Test
  using LinearAlgebra
  using IterativeSolvers
  using FillArrays


  function generate_stiffness_matrices(mh, fespaces, qdegree, α)
    tests,trials=fespaces
    Gridap.Helpers.@check num_levels(mh)==length(tests)
    Gridap.Helpers.@check num_levels(mh)==length(trials)
    matrices=Vector{PSparseMatrix}(undef,num_levels(mh))
    for i=1:num_levels(mh)
      model = get_level_model(mh,i)
      if (GridapP4est.i_am_in(model.parts))
        Ω  = Triangulation(model.dmodel)
        dΩ = Measure(Ω,qdegree)
        Vh = get_level_fe_space(tests[i])
        Uh = get_level_fe_space(trials[i])
        a(u,v)=∫(v⋅u)dΩ+∫((α*divergence(v))*divergence(u))dΩ
        A = assemble_matrix(a,Uh,Vh)
        # map_parts(A.owned_owned_values) do Ah
        #   println(eigvals(Array(Ah)))
        # end
        matrices[i]=A
      end
    end
    matrices
  end

  function generate_patch_based_smoothers(mh,fespaces,order,α)
    tests,trials=fespaces
    Gridap.Helpers.@check num_levels(mh)==length(tests)
    Gridap.Helpers.@check num_levels(mh)==length(trials)
    smoothers=Vector{RichardsonSmoother}(undef,num_levels(mh)-1)
    reffe=ReferenceFE(raviart_thomas,Float64,order)
    for i=1:num_levels(mh)-1
      model = get_level_model(mh,i)
      if (GridapP4est.i_am_in(model.parts))
        Vh = get_level_fe_space(tests[i])
        Uh = get_level_fe_space(trials[i])
        PD=PatchDecomposition(model.dmodel)
        #map_parts(PD.patch_decompositions,Vh.spaces) do PD, Vh
        #  println(PD.patch_cells)
        #  println(PD.patch_cells_faces_on_boundary[2])
        #  println(get_cell_dof_ids(Vh))
        #end
        Ph=PatchFESpace(model.dmodel,reffe,DivConformity(),PD,Vh)
        #map_parts(Ph.spaces) do space
        #  println(get_cell_dof_ids(space))
        #end
        Ω  = Triangulation(PD)
        dΩ = Measure(Ω,2*(order+1))
        a(u,v)=∫(v⋅u)dΩ+∫((α*divergence(v))*divergence(u))dΩ
        PLS=PatchBasedLinearSolver(a,Ph,LUSolver())
        smoothers[i]=RichardsonSmoother(PLS,1,1.0/3.0)
      end
    end
    smoothers
  end

  u(x) = VectorValue(x[1],x[2])
  f(x) = VectorValue(2.0*x[2]*(1.0-x[1]*x[1]),2.0*x[1]*(1-x[2]*x[2]))

  function run(parts,
               coarse_grid_partition,
               num_parts_x_level,
               order,
               α)
    domain=(0,1,0,1)

    reffe=ReferenceFE(raviart_thomas,Float64,order)
    qdegree=2*(order+1)

    cmodel=CartesianDiscreteModel(domain,coarse_grid_partition)
    mh=ModelHierarchy(parts,cmodel,num_parts_x_level)

    tests    = TestFESpace(mh,reffe)
    trials   = TrialFESpace(u,tests)
    fespaces = (tests,trials)
    smatrices= generate_stiffness_matrices(mh,fespaces,qdegree,α)
    smoothers= generate_patch_based_smoothers(mh,fespaces,order,α)
    interp,restrict=setup_interpolations_and_restrictions(mh,fespaces,qdegree)

    model=get_level_model(mh,1)
    Ω = Triangulation(model.dmodel)
    dΩ = Measure(Ω,qdegree)
    Vh = tests[1]
    Uh = trials[1]
    a(u,v)=∫(v⋅u)dΩ+∫((α*divergence(v))*divergence(u))dΩ
    l(v)=∫(v⋅f)dΩ
    Vh = get_level_fe_space(tests[1])
    Uh = get_level_fe_space(trials[1])

    op=AffineFEOperator(a,l,Uh,Vh)
    A=op.op.matrix
    b=op.op.vector
    x=PVector(0.0,A.cols)
    solver=GMGLinearSolver(mh,
                    smatrices,
                    interp,
                    restrict,
                    pre_smoothers=smoothers,
                    post_smoothers=smoothers,
                    maxiter=1,
                    rtol=1.0e-06,
                    verbose=false,
                    mode=:preconditioner)
    ss=symbolic_setup(solver,A)
    ns=numerical_setup(ss,A)

    x,history=IterativeSolvers.cg!(x,
                         A,
                         b;
                         verbose=i_am_main(parts),
                         reltol=1.0e-06,
                         Pl=ns,
                         log=true)

    # uh=FEFunction(Uh,x)
    # Error norms and print solution
    # e = u-uh
    # e_l2 = sum(∫(e*e)dΩ)
    # tol = 1.0e-9
    # @test e_l2 < tol
    # map_parts(parts) do part
    #  if (part==1)
    #    println("$(e_l2) < $(tol)\n")
    #  end
    # end
    model_hierarchy_free!(mh)
    history.iters,num_free_dofs(Vh)
  end
  if !MPI.Initialized()
    MPI.Init()
  end
  parts = get_part_ids(mpi,1)
  order=0

  num_refinements=[1,2,3,4,5]
  alpha_exps=[0,1,2,3,4]
  iter_matrix=zeros(Int,5,5)
  coarse_grid_partition=(5,5)
  free_dofs=Vector{Int64}(undef,length(num_refinements))

  for ref=1:length(num_refinements)
      num_parts_x_level=[1 for i=1:num_refinements[ref]+1]
      for alpha_exp=1:length(alpha_exps)
        α=10.0^alpha_exps[alpha_exp]
        num_iters,num_free_dofs=run(parts,coarse_grid_partition,num_parts_x_level,order,α)
        free_dofs[ref]=num_free_dofs
        iter_matrix[ref,alpha_exp]=num_iters
      end
  end
  println(iter_matrix)
  println(free_dofs)

  MPI.Finalize()
end
