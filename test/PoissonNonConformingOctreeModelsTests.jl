module PoissonNonConformingOctreeModelsTests
  using P4est_wrapper
  using GridapP4est
  using Gridap
  using PartitionedArrays
  using GridapDistributed
  using MPI
  using Gridap.FESpaces
  using FillArrays
  using Logging

  include("CoarseDiscreteModelsTools.jl")

  function generate_analytical_problem_functions(T::Type{Float64},order)
    # Define manufactured functions
    u(x) = x[1]+x[2]^order
    f(x) = -Δ(u)(x)
    u,f
  end

  function generate_analytical_problem_functions(T::Type{VectorValue{2,Float64}},order)
    # Define manufactured functions
    u(x) = VectorValue(x[1]+x[2]^order,x[1]^order+x[2])
    f(x) = -Δ(u)(x)
    u,f
  end
  
  function generate_analytical_problem_functions(T::Type{VectorValue{3,Float64}},order)
    # Define manufactured functions
    u(x) = VectorValue(x[1]+x[2]^order+x[3],x[1]^order+x[2]+x[3],x[1]+x[2]+x[3]^order)
    f(x) = -Δ(u)(x)
    u,f
  end

  function test_transfer_ops_and_redistribute(ranks,dmodel,order,cg_or_dg,T::Type)
    @assert cg_or_dg == :cg || cg_or_dg == :dg 
    if (cg_or_dg==:dg)
      @assert T==Float64
    end 

    conformity=cg_or_dg==:cg ? :H1 : :L2
    
    # Define manufactured functions
    u,f = generate_analytical_problem_functions(T,order)
    degree = 2*order+1
    reffe=ReferenceFE(lagrangian,T,order)
    if (cg_or_dg == :cg)
      VH=FESpace(dmodel,reffe,conformity=conformity;dirichlet_tags="boundary")
    else 
      VH=FESpace(dmodel,reffe,conformity=conformity)
    end 
    UH=TrialFESpace(VH,u)
    ref_coarse_flags=map(ranks,partition(get_cell_gids(dmodel.dmodel))) do rank,indices
        flags=zeros(Cint,length(indices))
        flags.=nothing_flag
        
        flags[1]=refine_flag
        flags[own_length(indices)]=refine_flag

        # To create some unbalance
        if (rank%2==0 && own_length(indices)>1)
            flags[div(own_length(indices),2)]=refine_flag
        end
        flags
    end 
    fmodel,glue=Gridap.Adaptivity.adapt(dmodel,ref_coarse_flags);
    # map(ranks,glue) do rank, glue 
    #   if rank==2
    #     print(glue.n2o_faces_map[end]); print("\n")
    #   end  
    # end
    if (cg_or_dg == :cg)
      Vh=FESpace(fmodel,reffe,conformity=conformity;dirichlet_tags="boundary")
    else
      Vh=FESpace(fmodel,reffe,conformity=conformity)
    end 
    Uh=TrialFESpace(Vh,u)

    ΩH  = Triangulation(dmodel)
    dΩH = Measure(ΩH,degree)

    if (cg_or_dg==:cg)
      aHcg(u,v) = ∫( ∇(v)⊙∇(u) )*dΩH
      bHcg(v) = ∫(v⋅f)*dΩH
      op = AffineFEOperator(aHcg,bHcg,UH,VH)
    else
      h = 2
      γ = 10

      ΛH = Skeleton(dmodel)
      ΓH = Boundary(dmodel,tags="boundary")
      
      dΛH = Measure(ΛH,2*order)
      dΓH = Measure(ΓH,2*order)
      
      n_ΓH = get_normal_vector(ΓH)
      n_ΛH = get_normal_vector(ΛH)

      aHdg(u,v) =
          ∫( ∇(v)⋅∇(u) )*dΩH +
          ∫( (γ/h)*v*u  - v*(n_ΓH⋅∇(u)) - (n_ΓH⋅∇(v))*u )*dΓH +
          ∫( (γ/h)*jump(v*n_ΛH)⋅jump(u*n_ΛH) -
          jump(v*n_ΛH)⋅mean(∇(u)) -
          mean(∇(v))⋅jump(u*n_ΛH) )*dΛH

      bHdg(v) =
          ∫( v*f )*dΩH +
          ∫( (γ/h)*v*u - (n_ΓH⋅∇(v))*u )*dΓH
      op = AffineFEOperator(aHdg,bHdg,UH,VH)
    end 

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

    if (cg_or_dg==:cg)
      ahcg(u,v) = ∫( ∇(v)⊙∇(u) )*dΩh
      bhcg(v) = ∫(v⋅f)*dΩh
      op = AffineFEOperator(ahcg,bhcg,Uh,Vh)
    else
      h = 2
      γ = 10

      Λh = Skeleton(fmodel)
      Γh = Boundary(fmodel,tags="boundary")
      
      dΛh = Measure(Λh,2*order)
      dΓh = Measure(Γh,2*order)
      
      n_Γh = get_normal_vector(Γh)
      n_Λh = get_normal_vector(Λh)

      ahdg(u,v) =
          ∫( ∇(v)⋅∇(u) )*dΩh +
          ∫( (γ/h)*v*u  - v*(n_Γh⋅∇(u)) - (n_Γh⋅∇(v))*u )*dΓh +
          ∫( (γ/h)*jump(v*n_Λh)⋅jump(u*n_Λh) -
          jump(v*n_Λh)⋅mean(∇(u)) -
          mean(∇(v))⋅jump(u*n_Λh) )*dΛh

      bhdg(v) =
          ∫( v*f )*dΩh +
          ∫( (γ/h)*v*u - (n_Γh⋅∇(v))*u )*dΓh
      op = AffineFEOperator(ahdg,bhdg,Uh,Vh)
    end 

    uh = solve(op)
    e = u - uh

    writevtk(ΩH, "ctrian", cellfields=["uH"=>uH])
    writevtk(Ωh, "ftrian", cellfields=["uh"=>uh])

    # # Compute errors

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

    weights=map(ranks,fmodel.dmodel.models) do rank,lmodel
      if (rank%2==0)
        zeros(Cint,num_cells(lmodel))
      else
        ones(Cint,num_cells(lmodel))
      end
    end 
    fmodel_red, red_glue=GridapDistributed.redistribute(fmodel);
    if (cg_or_dg==:cg)
      Vhred=FESpace(fmodel_red,reffe,conformity=conformity;dirichlet_tags="boundary")
      Uhred=TrialFESpace(Vhred,u)
    else
      Vhred=FESpace(fmodel_red,reffe,conformity=conformity)
      Uhred=TrialFESpace(Vhred,u)
    end 

    Ωhred  = Triangulation(fmodel_red)
    dΩhred = Measure(Ωhred,degree)

    if (cg_or_dg==:cg)
      ahcgred(u,v) = ∫( ∇(v)⊙∇(u) )*dΩhred
      bhcgred(v)   = ∫(v⋅f)*dΩhred
      op = AffineFEOperator(ahcgred,bhcgred,Uhred,Vhred)
    else
      h = 2
      γ = 10

      Λhred = Skeleton(fmodel_red)
      Γhred = Boundary(fmodel_red,tags="boundary")
      
      dΛhred = Measure(Λhred,2*order)
      dΓhred = Measure(Γhred,2*order)
      
      n_Γhred = get_normal_vector(Γhred)
      n_Λhred = get_normal_vector(Λhred)

      ahdgred(u,v) =
          ∫( ∇(v)⋅∇(u) )*dΩhred +
          ∫( (γ/h)*v*u  - v*(n_Γhred⋅∇(u)) - (n_Γhred⋅∇(v))*u )*dΓhred +
          ∫( (γ/h)*jump(v*n_Λhred)⋅jump(u*n_Λhred) -
          jump(v*n_Λhred)⋅mean(∇(u)) -
          mean(∇(v))⋅jump(u*n_Λhred) )*dΛhred

      bhdgred(v) =
          ∫( v*f )*dΩhred +
          ∫( (γ/h)*v*u - (n_Γhred⋅∇(v))*u )*dΓhred
      op = AffineFEOperator(ahdgred,bhdgred,Uhred,Vhred)
    end

    
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
            cg_or_dg,
            T::Type) where Dc

    @assert cg_or_dg == :cg || cg_or_dg == :dg 
    if (cg_or_dg==:dg)
      @assert T==Float64
    end 

    # Define manufactured functions
    u,f = generate_analytical_problem_functions(T,order)

    degree = 2*order+1
    ref_coarse_flags=map(ranks,partition(get_cell_gids(dmodel.dmodel))) do rank,indices
        flags=zeros(Int,length(indices))
        flags.=nothing_flag        
        if (rank==1)
           flags[1:2^Dc].=coarsen_flag
        else 
           flags[own_length(indices)]=refine_flag
        end 
        flags
    end
    fmodel,glue=Gridap.Adaptivity.adapt(dmodel,ref_coarse_flags);

    conformity=cg_or_dg==:cg ? :H1 : :L2

    reffe=ReferenceFE(lagrangian,T,order)
    if (conformity==:L2)
      VH=FESpace(dmodel,reffe,conformity=conformity)
    else 
      VH=FESpace(dmodel,reffe,conformity=conformity;dirichlet_tags="boundary")
    end 
    UH=TrialFESpace(VH,u)

    if (conformity==:L2)
      Vh=FESpace(fmodel,reffe,conformity=conformity)
    else   
      Vh=FESpace(fmodel,reffe,conformity=conformity;dirichlet_tags="boundary")
    end 
    Uh=TrialFESpace(Vh,u)
    ΩH  = Triangulation(dmodel)
    dΩH = Measure(ΩH,degree)

    if (cg_or_dg==:cg)
      aHcg(u,v) = ∫( ∇(v)⊙∇(u) )*dΩH
      bHcg(v) = ∫(v⋅f)*dΩH
      op = AffineFEOperator(aHcg,bHcg,UH,VH)
    else
      h = 2
      γ = 10

      ΛH = Skeleton(dmodel)
      ΓH = Boundary(dmodel,tags="boundary")
      
      dΛH = Measure(ΛH,2*order)
      dΓH = Measure(ΓH,2*order)
      
      n_ΓH = get_normal_vector(ΓH)
      n_ΛH = get_normal_vector(ΛH)

      aHdg(u,v) =
         ∫( ∇(v)⋅∇(u) )*dΩH +
         ∫( (γ/h)*v*u  - v*(n_ΓH⋅∇(u)) - (n_ΓH⋅∇(v))*u )*dΓH +
         ∫( (γ/h)*jump(v*n_ΛH)⋅jump(u*n_ΛH) -
         jump(v*n_ΛH)⋅mean(∇(u)) -
         mean(∇(v))⋅jump(u*n_ΛH) )*dΛH

      bHdg(v) =
        ∫( v*f )*dΩH +
        ∫( (γ/h)*v*u - (n_ΓH⋅∇(v))*u )*dΓH
        op = AffineFEOperator(aHdg,bHdg,UH,VH)
    end 

    uH = solve(op)
    e = u - uH

    # # Compute errors
    el2 = sqrt(sum( ∫( e⋅e )*dΩH ))
    eh1 = sqrt(sum( ∫( e⋅e + ∇(e)⊙∇(e) )*dΩH ))

    tol=1e-5
    @assert el2 < tol
    @assert eh1 < tol

    Ωh  = Triangulation(fmodel)
    dΩh = Measure(Ωh,degree)

    if (cg_or_dg==:cg)
      ahcg(u,v) = ∫( ∇(v)⊙∇(u) )*dΩh
      bhcg(v) = ∫(v⋅f)*dΩh
      op = AffineFEOperator(ahcg,bhcg,Uh,Vh)
    else
      h = 2
      γ = 10

      Λh = Skeleton(fmodel)
      Γh = Boundary(fmodel,tags="boundary")
      
      dΛh = Measure(Λh,2*order)
      dΓh = Measure(Γh,2*order)
      
      n_Γh = get_normal_vector(Γh)
      n_Λh = get_normal_vector(Λh)

      ahdg(u,v) =
          ∫( ∇(v)⋅∇(u) )*dΩh +
          ∫( (γ/h)*v*u  - v*(n_Γh⋅∇(u)) - (n_Γh⋅∇(v))*u )*dΓh +
          # with p=4 MPI tasks, (γ/h)*jump(v*n_Λh)⋅jump(u*n_Λh) fails!
          # did not have time you to understand the cause, workaround 
          # is the line below
          ∫( (γ/h)*(jump(v*n_Λh)⋅jump(u*n_Λh)) -
             jump(v*n_Λh)⋅mean(∇(u)) -
             mean(∇(v))⋅jump(u*n_Λh) )*dΛh

      bhdg(v) =
          ∫( v*f )*dΩh +
          ∫( (γ/h)*v*u - (n_Γh⋅∇(v))*u )*dΓh
      op = AffineFEOperator(ahdg,bhdg,Uh,Vh)
    end 

    uh = solve(op)
    e = u - uh

    # # Compute errors
    el2 = sqrt(sum( ∫( e⋅e )*dΩh ))
    eh1 = sqrt(sum( ∫( e⋅e + ∇(e)⊙∇(e) )*dΩh ))

    @assert el2 < tol
    @assert eh1 < tol

    # prolongation via interpolation
    uHh=interpolate(uH,Uh)
    e = uh - uHh
    el2 = sqrt(sum( ∫( e⋅e )*dΩh ))
    @assert el2 < tol
  end

  function test_2d(ranks,order,cg_or_dg,T::Type;num_amr_steps=5)
    coarse_model=CartesianDiscreteModel((0,1,0,1),(1,1))
    dmodel=OctreeDistributedDiscreteModel(ranks,coarse_model,2)
    test_refine_and_coarsen_at_once(ranks,dmodel,order,cg_or_dg,T)
    rdmodel=dmodel
    for i=1:num_amr_steps
     rdmodel=test_transfer_ops_and_redistribute(ranks,rdmodel,order,cg_or_dg,T)
    end
  end 

  function test_3d(ranks,order,cg_or_dg,T::Type;num_amr_steps=5)
    coarse_model=CartesianDiscreteModel((0,1,0,1,0,1),(1,1,1))
    dmodel=OctreeDistributedDiscreteModel(ranks,coarse_model,2)
    test_refine_and_coarsen_at_once(ranks,dmodel,order,cg_or_dg,T)
    rdmodel=dmodel
    for i=1:num_amr_steps
      rdmodel=test_transfer_ops_and_redistribute(ranks,rdmodel,order,cg_or_dg,T)
    end
  end 

  function test(ranks,TVDc::Type{Val{Dc}}, perm, order,cg_or_dg,T::Type) where Dc
    coarse_model = setup_model(TVDc,perm)
    model = OctreeDistributedDiscreteModel(ranks, coarse_model, 1)
    test_transfer_ops_and_redistribute(ranks,model,order,cg_or_dg,T)
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
    #debug_logger = ConsoleLogger(stderr, Logging.Debug)
    #global_logger(debug_logger); # Enable the debug logger globally
    ranks = distribute(LinearIndices((MPI.Comm_size(MPI.COMM_WORLD),)))
    for Dc=2:3, perm in (1,2,4), order=(1,2), scalar_or_vector in (:scalar,)
      test(ranks,Val{Dc},perm,order,:dg,_field_type(Val{Dc}(),scalar_or_vector))
    end
    # for Dc=2:3, perm in (1,2), order in (1,4), scalar_or_vector in (:vector,)
    #  test(ranks,Val{Dc},perm,order,:cg,_field_type(Val{Dc}(),scalar_or_vector))
    # end
    for order=2:2, scalar_or_vector in (:scalar,)
      test_2d(ranks,order,:dg,_field_type(Val{2}(),scalar_or_vector), num_amr_steps=5)
      test_3d(ranks,order,:dg,_field_type(Val{3}(),scalar_or_vector), num_amr_steps=4)
    end
    # for order=2:2, scalar_or_vector in (:scalar,:vector)
    #  test_2d(ranks,order,:cg,_field_type(Val{2}(),scalar_or_vector), num_amr_steps=5)
    #  test_3d(ranks,order,:cg,_field_type(Val{3}(),scalar_or_vector), num_amr_steps=4)
    # end
  end
end
