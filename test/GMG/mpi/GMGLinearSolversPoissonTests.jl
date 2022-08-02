module GMGLinearSolverPoissonTests
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


  function generate_stiffness_matrices(mh, fespaces, qdegree)
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
        a(u,v)=∫(∇(v)⋅∇(u))dΩ
        A = assemble_matrix(a,Uh,Vh)
        matrices[i]=A
      end
    end
    matrices
  end

  function generate_patch_based_smoothers(mh,fespaces,order)
    tests,trials=fespaces
    Gridap.Helpers.@check num_levels(mh)==length(tests)
    Gridap.Helpers.@check num_levels(mh)==length(trials)
    smoothers=Vector{RichardsonSmoother}(undef,num_levels(mh)-1)
    reffe=ReferenceFE(lagrangian,Float64,order)
    for i=1:num_levels(mh)-1
      model = get_level_model(mh,i)
      if (GridapP4est.i_am_in(model.parts))
        Vh = get_level_fe_space(tests[i])
        Uh = get_level_fe_space(trials[i])
        PD=PatchDecomposition(model.dmodel)
        Ph=PatchFESpace(model.dmodel,reffe,H1Conformity(),PD,Vh)
        Ω  = Triangulation(PD)
        dΩ = Measure(Ω,2*order+1)
        a(u,v)=∫(∇(v)⋅∇(u))dΩ
        PLS=PatchBasedLinearSolver(a,Ph,LUSolver())
        smoothers[i]=RichardsonSmoother(PLS,10)
      end
    end
    smoothers
  end

  # Manufactured solution
  u(x) = x[1] + x[2]
  f(x) = -Δ(u)(x)

  function run(parts,
               coarse_grid_partition,
               num_parts_x_level,
               order)
    domain=(0,1,0,1)

    reffe=ReferenceFE(lagrangian,Float64,order)
    qdegree=2*order+1

    cmodel=CartesianDiscreteModel(domain,coarse_grid_partition)
    mh=ModelHierarchy(parts,cmodel,num_parts_x_level)

    tests    = TestFESpace(mh,reffe,dirichlet_tags="boundary")
    trials   = TrialFESpace(u,tests)
    fespaces = (tests,trials)
    smatrices= generate_stiffness_matrices(mh,fespaces,qdegree)
    # smoothers= generate_patch_based_smoothers(mh,fespaces,order)
    smoothers = Fill(RichardsonSmoother(JacobiLinearSolver(),10),num_levels(mh)-1)
    interp,restrict=setup_interpolations_and_restrictions(mh,fespaces,qdegree)

    model=get_level_model(mh,1)
    Ω = Triangulation(model.dmodel)
    dΩ = Measure(Ω,qdegree)
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

    IterativeSolvers.cg!(x,
                         A,
                         b;
                         verbose=i_am_main(parts),
                         reltol=1.0e-06,
                         Pl=ns)

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
  parts = get_part_ids(mpi,4)
  order=1
  num_parts_x_level=[4,4,2,1]
  coarse_grid_partition=(2,2)
  run(parts,coarse_grid_partition,num_parts_x_level,order)
  MPI.Finalize()
end
