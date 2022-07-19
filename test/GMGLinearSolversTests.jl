module GMGLinearSolverTests
  using MPI
  using Gridap
  using PartitionedArrays
  using GridapDistributed
  using GridapP4est
  using P4est_wrapper
  using Test
  using LinearAlgebra
  using IterativeSolvers


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

  function setup_interpolations_and_restrictions(mh,fespaces,qdegree)
    nlevs = num_levels(mh)
    interpolations=Vector{InterpolationMat}(undef,nlevs-1)
    restrictions=Vector{RestrictionMat}(undef,nlevs-1)
    for l=1:nlevs-1
      model = get_level_model(mh,l)
      if (GridapP4est.i_am_in(model.parts))
        interpolations[l]=InterpolationMat(mh,fespaces,l,qdegree)
        restrictions[l]=RestrictionMat(mh,fespaces,l,qdegree)
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

    order=1
    reffe=ReferenceFE(lagrangian,Float64,order)
    qdegree=2*order+1

    num_parts_x_level = [4,4,2,1]
    cmodel=CartesianDiscreteModel(domain,(2,2))
    mh=ModelHierarchy(parts,cmodel,num_parts_x_level)

    tests    = TestFESpace(mh,reffe,dirichlet_tags="boundary")
    trials   = TrialFESpace(u,tests)
    fespaces = (tests,trials)
    smatrices= generate_stiffness_matrices(mh,fespaces,qdegree)
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



    GMG!(x,
         b,
         mh,
         smatrices,
         interp,
         restrict;
         rtol=1.0e-06,
         maxiter=200,
         smoother=JacobiSmoother(10))

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
  run(parts,(1,1))
  MPI.Finalize()
end
