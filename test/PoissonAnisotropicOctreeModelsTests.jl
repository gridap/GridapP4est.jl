module PoissonAnisotropicOctreeModelsTests
  using P4est_wrapper
  using GridapP4est
  using Gridap
  using PartitionedArrays
  using GridapDistributed
  using MPI
  using Gridap.FESpaces
  using FillArrays
  using Logging
  using LinearAlgebra

  include("CoarseDiscreteModelsTools.jl")

  function generate_analytical_problem_functions(T::Type{Float64},order)
    # Define manufactured functions
    u(x) = x[1]+x[2]^order
    f(x) = -Δ(u)(x)
    u,f
  end
  
  function generate_analytical_problem_functions(T::Type{VectorValue{3,Float64}},order)
    # Define manufactured functions
    u(x) = VectorValue(x[1]+x[2]^order+x[3],x[1]^order+x[2]+x[3],x[1]+x[2]+x[3]^order)
    f(x) = -Δ(u)(x)
    u,f
  end

  function test_transfer_ops_and_redistribute(ranks,dmodel,order,T::Type)
    # Define manufactured functions
    u,f = generate_analytical_problem_functions(T,order)
    degree = 2*order+1
    reffe=ReferenceFE(lagrangian,T,order)
    VH=FESpace(dmodel,reffe;dirichlet_tags="boundary")
    UH=TrialFESpace(VH,u)
    num_local_cols=GridapP4est.num_locally_owned_columns(dmodel)
    ref_coarse_flags=map(ranks,num_local_cols) do rank,num_local_cols
      println("[rank:$(rank)] $(num_local_cols)")
      
      flags=zeros(Cint,num_local_cols)
      flags.=nothing_flag
      
      flags[1]=refine_flag
      flags[end]=refine_flag

      # To create some unbalance
      if (rank%2==0)
         flags[div(num_local_cols,2)]=refine_flag
      end
      flags
    end
    fmodel,glue=GridapP4est.horizontally_adapt(dmodel,ref_coarse_flags);

    writevtk(fmodel,"fmodel")

    Vh=FESpace(fmodel,reffe,conformity=:H1;dirichlet_tags="boundary")
    Uh=TrialFESpace(Vh,u)

    ΩH  = Triangulation(dmodel)
    dΩH = Measure(ΩH,degree)

    aH(u,v) = ∫( ∇(v)⊙∇(u) )*dΩH
    bH(v) = ∫(v⋅f)*dΩH

    op = AffineFEOperator(aH,bH,UH,VH)
    uH = solve(op)
    e = u - uH

    # # Compute errors
    el2 = sqrt(sum( ∫( e⋅e )*dΩH ))
    eh1 = sqrt(sum( ∫( e⋅e + ∇(e)⊙∇(e) )*dΩH ))

    tol=1e-5
    println("[SOLVE COARSE] el2 < tol: $(el2) < $(tol)")
    println("[SOLVE COARSE] eh1 < tol: $(eh1) < $(tol)")
    @assert el2 < tol
    @assert eh1 < tol

    Ωh  = Triangulation(fmodel)
    dΩh = Measure(Ωh,degree)

    ah(u,v) = ∫( ∇(v)⊙∇(u) )*dΩh
    bh(v) = ∫(v⋅f)*dΩh

    op  = AffineFEOperator(ah,bh,Uh,Vh)
    uh = solve(op)

    uh2dofs=get_free_dof_values(interpolate(u, Uh))

    # println(op.op.matrix.matrix_partition.item_ref[]*uh2dofs.vector_partition.item_ref[]-
    #         op.op.vector.vector_partition.item_ref[])
    # writevtk(ΩH, "ctrian", cellfields=["uH"=>uH])
    # writevtk(Ωh, "ftrian", cellfields=["uh"=>uh])

    # # Compute errors
    e = u - uh
    el2 = sqrt(sum( ∫( e⋅e )*dΩh ))
    eh1 = sqrt(sum( ∫( e⋅e + ∇(e)⊙∇(e) )*dΩh )) 
    println("[SOLVE FINE] el2 < tol: $(el2) < $(tol)")
    println("[SOLVE FINE] eh1 < tol: $(eh1) < $(tol)")
    @assert el2 < tol
    @assert eh1 < tol

    # prolongation via interpolation
    uHh=interpolate(uH,Uh)   
    e = uh - uHh
    el2 = sqrt(sum( ∫( e⋅e )*dΩh ))
    println("[INTERPOLATION] el2 < tol: $(el2) < $(tol)")
    @assert el2 < tol

    # prolongation via L2-projection 
    # Coarse FEFunction -> Fine FEFunction, by projection
    ahp(u,v)  = ∫(v⋅u)*dΩh
    lhp(v)    = ∫(v⋅uH)*dΩh
    oph      = AffineFEOperator(ahp,lhp,Uh,Vh)
    uHh      = solve(oph)
    e = uh - uHh
    el2 = sqrt(sum( ∫( e⋅e )*dΩh ))
    println("[L2 PROJECTION] el2 < tol: $(el2) < $(tol)")
    @assert el2 < tol

    # restriction via interpolation
    uhH=interpolate(uh,UH) 
    e = uH - uhH
    el2 = sqrt(sum( ∫( e⋅e )*dΩh ))
    println("[INTERPOLATION] el2 < tol: $(el2) < $(tol)")
    @assert el2 < tol

    # restriction via L2-projection
    dΩhH = Measure(ΩH,Ωh,2*order)
    aHp(u,v) = ∫(v⋅u)*dΩH
    lHp(v)   = ∫(v⋅uh)*dΩhH
    oph     = AffineFEOperator(aHp,lHp,UH,VH)
    uhH     = solve(oph)
    e       = uH - uhH
    el2     = sqrt(sum( ∫( e⋅e )*dΩH ))

    fmodel_red, red_glue=GridapDistributed.redistribute(fmodel);
    Vhred=FESpace(fmodel_red,reffe,conformity=:H1;dirichlet_tags="boundary")
    Uhred=TrialFESpace(Vhred,u)

    Ωhred  = Triangulation(fmodel_red)
    dΩhred = Measure(Ωhred,degree)

    ahred(u,v) = ∫( ∇(v)⊙∇(u) )*dΩhred
    bhred(v)   = ∫(v⋅f)*dΩhred

    op    = AffineFEOperator(ahred,bhred,Uhred,Vhred)
    uhred = solve(op)
    e = u - uhred
    el2 = sqrt(sum( ∫( e⋅e )*dΩhred ))
    println("[SOLVE FINE REDISTRIBUTED] el2 < tol: $(el2) < $(tol)")
    @assert el2 < tol

    uhred2 = GridapDistributed.redistribute_fe_function(uh,Vhred,fmodel_red,red_glue)
    e = u - uhred2

    el2 = sqrt(sum( ∫( e⋅e )*dΩhred ))
    println("[REDISTRIBUTE SOLUTION] el2 < tol: $(el2) < $(tol)")
    @assert el2 < tol

     fmodel_red
  end

  function test_refine_and_coarsen_at_once(ranks,
            dmodel::OctreeDistributedDiscreteModel{Dc},
            order,
            T::Type) where Dc

    # Define manufactured functions
    u,f = generate_analytical_problem_functions(T,order)

    num_local_cols=GridapP4est.num_locally_owned_columns(dmodel)

    degree = 2*order+1
    ref_coarse_flags=map(ranks,num_local_cols) do rank,num_local_cols
        flags=zeros(Int,num_local_cols)
        flags.=nothing_flag        
        #if (rank==1)
        #   flags[1:4].=coarsen_flag
        #else 
        flags[1]=refine_flag
        #end 
        flags
    end
    fmodel,glue=GridapP4est.horizontally_adapt(dmodel,ref_coarse_flags);

    reffe=ReferenceFE(lagrangian,T,order)
    VH=FESpace(dmodel,reffe,conformity=:H1;dirichlet_tags="boundary")
    UH=TrialFESpace(VH,u)

    Vh=FESpace(fmodel,reffe,conformity=:H1;dirichlet_tags="boundary")
    Uh=TrialFESpace(Vh,u)
    ΩH  = Triangulation(dmodel)
    dΩH = Measure(ΩH,degree)

    aH(u,v) = ∫( ∇(v)⊙∇(u) )*dΩH
    bH(v) = ∫(v⋅f)*dΩH

    op = AffineFEOperator(aH,bH,UH,VH)
    uH = solve(op)
    e = u - uH

    # # Compute errors
    el2 = sqrt(sum( ∫( e⋅e )*dΩH ))
    eh1 = sqrt(sum( ∫( e⋅e + ∇(e)⊙∇(e) )*dΩH ))

    tol=1e-5
    println("[SOLVE INITIAL MESH] el2 < tol: $(el2) < $(tol)")
    println("[SOLVE INITIAL MESH] eh1 < tol: $(eh1) < $(tol)")
    @assert el2 < tol
    @assert eh1 < tol

    Ωh  = Triangulation(fmodel)
    dΩh = Measure(Ωh,degree)

    ah(u,v) = ∫( ∇(v)⊙∇(u) )*dΩh
    bh(v) = ∫(v⋅f)*dΩh

    op = AffineFEOperator(ah,bh,Uh,Vh)
    uh = solve(op)
    e = u - uh

    # # Compute errors
    el2 = sqrt(sum( ∫( e⋅e )*dΩh ))
    eh1 = sqrt(sum( ∫( e⋅e + ∇(e)⊙∇(e) )*dΩh ))
    println("[SOLVE ADAPTED MESH] el2 < tol: $(el2) < $(tol)")
    println("[SOLVE ADAPTED MESH] eh1 < tol: $(eh1) < $(tol)")
    @assert el2 < tol
    @assert eh1 < tol

    # prolongation via interpolation
    uHh=interpolate(uH,Uh)
    e = uh - uHh
    el2 = sqrt(sum( ∫( e⋅e )*dΩh ))
    println("[INTERPOLATION INITIAL-ADAPTED MESH] el2 < tol: $(el2) < $(tol)")
    @assert el2 < tol
  end

  function test_3d(ranks,order,T::Type;num_amr_steps=5)
    coarse_model=CartesianDiscreteModel((0,1,0,1),(1,1))
    dmodel=AnisotropicallyAdapted3DDistributedDiscreteModel(ranks, coarse_model, 2, 2)
    #writevtk(dmodel.dmodel,"model")
    #test_refine_and_coarsen_at_once(ranks,dmodel,order,T)
    rdmodel=dmodel
    for i=1:num_amr_steps
     rdmodel=test_transfer_ops_and_redistribute(ranks,rdmodel,order,T)
     if i==1
      rdmodel,_=vertically_uniformly_refine(dmodel)
     end
    end
  end 

  function test(ranks,perm, order,T::Type)
    coarse_model = setup_model(Val{2},perm)
    model = AnisotropicallyAdapted3DDistributedDiscreteModel(ranks, coarse_model, 1, 1)
    test_transfer_ops_and_redistribute(ranks,model,order,T)
  end
  
  function _field_type(::Val{Dc}, scalar_or_vector::Symbol) where Dc
    if scalar_or_vector==:scalar
      Float64
    else 
      @assert scalar_or_vector==:vector
      VectorValue{Dc,Float64}
    end
  end 
  function run(distribute)
    # debug_logger = ConsoleLogger(stderr, Logging.Debug)
    # global_logger(debug_logger); # Enable the debug logger globally
    ranks = distribute(LinearIndices((MPI.Comm_size(MPI.COMM_WORLD),)))
    for perm=1:4, order=1:4, scalar_or_vector in (:scalar,)
      test(ranks,perm,order,_field_type(Val{3}(),scalar_or_vector))
    end
    for perm in (1,2), order in (1,4), scalar_or_vector in (:vector,)
      test(ranks,perm,order,_field_type(Val{3}(),scalar_or_vector))
    end
    for order=1:1, scalar_or_vector in (:scalar,:vector)
      test_3d(ranks,order,_field_type(Val{3}(),scalar_or_vector), num_amr_steps=4)
    end
  end
end

