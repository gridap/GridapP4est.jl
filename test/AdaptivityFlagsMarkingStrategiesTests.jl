module AdaptivityFlagsMarkingStrategiesTests
  using P4est_wrapper
  using GridapP4est
  using Gridap
  using PartitionedArrays
  using GridapDistributed
  using MPI
  using Gridap.FESpaces
  using FillArrays
  using Logging

  function test_refine_and_coarsen_at_once(ranks,
                                           dmodel::OctreeDistributedDiscreteModel{Dc},
                                           order) where Dc

    refinement_fraction=0.2
    coarsening_fraction=0.05
    cell_partition   = get_cell_gids(dmodel)

    degree = 2*order+1
    reffe=ReferenceFE(lagrangian,Float64,order)

    Vh=FESpace(dmodel,reffe,conformity=:H1)
    Uh=TrialFESpace(Vh)

    
    α  = 200
    r  = 0.7
    if (Dc==2)
      xc = VectorValue(-0.05, -0.05)
    else
      xc = VectorValue(-0.05, -0.05, -0.05)
    end 
    g(x) = atan(α*(sqrt((x-xc)·(x-xc))-r))

    gh=interpolate(g,Uh)
    Ω  = Triangulation(with_ghost,dmodel)
    dΩ = Measure(Ω,degree)
    dc=∫((g-gh)*(g-gh))dΩ
    error_indicators=map(local_views(dc)) do dc 
      sqrt.(get_array(dc))
    end

    ref_coarse_flags=map(ranks,partition(get_cell_gids(dmodel.dmodel))) do rank,indices
        flags=zeros(Cint,length(indices))
        flags.=nothing_flag        
    end

    adaptive_strategy=FixedFractionAdaptiveFlagsMarkingStrategy(refinement_fraction,coarsening_fraction)
    update_adaptivity_flags!(ref_coarse_flags,
                             adaptive_strategy,
                             partition(cell_partition),
                             error_indicators;
                             verbose=true)

    model,glue=adapt(dmodel,ref_coarse_flags);

    writevtk(model, "model$(Dc)")

    # Define manufactured functions
    u(x) = x[1]+x[2]^order
    f(x) = -Δ(u)(x)
    degree = 2*order+1
    reffe=ReferenceFE(lagrangian,Float64,order)

    Vh=FESpace(model,reffe,conformity=:H1;dirichlet_tags="boundary")
    Uh=TrialFESpace(Vh,u)

    Ω  = Triangulation(model)
    dΩ = Measure(Ω,degree)

    a(u,v) = ∫( ∇(v)⊙∇(u) )*dΩ
    b(v) = ∫(v*f)*dΩ

    op = AffineFEOperator(a,b,Uh,Vh)
    uh = solve(op)
    e = u - uh

    # # Compute errors
    el2 = sqrt(sum( ∫( e*e )*dΩ ))
    eh1 = sqrt(sum( ∫( e*e + ∇(e)⋅∇(e) )*dΩ ))

    tol=1e-6
    println("[SOLVE] el2 < tol: $(el2) < $(tol)")
    println("[SOLVE] eh1 < tol: $(eh1) < $(tol)")
    @assert el2 < tol
    @assert eh1 < tol

    model
  end

  function test_2d(ranks,order)
    coarse_model=CartesianDiscreteModel((0,1,0,1),(1,1))
    dmodel=OctreeDistributedDiscreteModel(ranks,coarse_model,2)
    dmodel=test_refine_and_coarsen_at_once(ranks,dmodel,order)
    # rdmodel=dmodel
    for i=1:4
      dmodel=test_refine_and_coarsen_at_once(ranks,dmodel,order)
    end
  end 

  function test_3d(ranks,order)
    coarse_model=CartesianDiscreteModel((0,1,0,1,0,1),(1,1,1))
    dmodel=OctreeDistributedDiscreteModel(ranks,coarse_model,2)
    dmodel=test_refine_and_coarsen_at_once(ranks,dmodel,order)
    for i=1:4
        dmodel=test_refine_and_coarsen_at_once(ranks,dmodel,order)
    end
  end

  function run(distribute)
    # debug_logger = ConsoleLogger(stderr, Logging.Debug)
    # global_logger(debug_logger); # Enable the debug logger globally

    ranks = distribute(LinearIndices((MPI.Comm_size(MPI.COMM_WORLD),)))
    for order=1:1
      test_2d(ranks,order)
      test_3d(ranks,order)
    end
  end

end