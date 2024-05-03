module MaxwellNonConformingOctreeModelsTests
  using P4est_wrapper
  using GridapP4est
  using Gridap
  using PartitionedArrays
  using GridapDistributed
  using MPI
  using Gridap.FESpaces
  using FillArrays
  using Logging
  using Test

  function test_transfer_ops_and_redistribute(ranks,
                                              dmodel::GridapDistributed.DistributedDiscreteModel{Dc},
                                              order) where Dc
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
    
    # Solve coarse
    uH,UH=solve_maxwell(dmodel,order)
    check_error_maxwell(dmodel,order,uH)

    # Solve fine
    uh,Uh=solve_maxwell(fmodel,order)
    check_error_maxwell(fmodel,order,uh)

    Ωh = Triangulation(fmodel)
    degree = 2*(order+1)
    dΩh = Measure(Ωh,degree)

    # prolongation via interpolation
    uHh=interpolate(uH,Uh)
    e = uh - uHh
    el2 = sqrt(sum( ∫( e⋅e )*dΩh ))
    tol=1e-6
    println("[INTERPOLATION] el2 < tol: $(el2) < $(tol)")
    @assert el2 < tol

    # prolongation via L2-projection 
    # Coarse FEFunction -> Fine FEFunction, by projection
    ahp(u,v)  = ∫(v⋅u)*dΩh
    lhp(v)    = ∫(v⋅uH)*dΩh
    oph      = AffineFEOperator(ahp,lhp,Uh,Uh)
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
    ΩH = Triangulation(dmodel)
    degree = 2*(order+1)
    dΩH = Measure(ΩH,degree)
    
    dΩhH = Measure(ΩH,Ωh,2*order)
    aHp(u,v) = ∫(v⋅u)*dΩH
    lHp(v)   = ∫(v⋅uh)*dΩhH
    oph     = AffineFEOperator(aHp,lHp,UH,UH)
    uhH     = solve(oph)
    e       = uH - uhH
    el2     = sqrt(sum( ∫( e⋅e )*dΩH ))

    fmodel_red, red_glue=GridapDistributed.redistribute(fmodel);
    uh_red,Uh_red=solve_maxwell(fmodel_red,order)
    check_error_maxwell(fmodel_red,order,uh_red)

    trian = Triangulation(fmodel_red)
    degree = 2*(order+1)
    dΩhred = Measure(trian,degree)

    u_ex, f_ex=get_analytical_functions(Dc)

    uhred2 = GridapDistributed.redistribute_fe_function(uh,Uh_red,fmodel_red,red_glue)
    
    e = u_ex - uhred2
    el2 = sqrt(sum( ∫( e⋅e )*dΩhred ))
    println("[REDISTRIBUTE SOLUTION] el2 < tol: $(el2) < $(tol)")
    @assert el2 < tol

    fmodel_red
  end

  function test_refine_and_coarsen_at_once(ranks,
            dmodel::OctreeDistributedDiscreteModel{Dc},
            order) where Dc
    degree = 2*order+1
    ref_coarse_flags=map(ranks,partition(get_cell_gids(dmodel.dmodel))) do rank,indices
        flags=zeros(Cint,length(indices))
        flags.=nothing_flag        
        if (rank==1)
          flags[1:min(2^Dc,own_length(indices))].=coarsen_flag
        end  
        flags[own_length(indices)]=refine_flag 
        flags
    end
    fmodel,glue=Gridap.Adaptivity.adapt(dmodel,ref_coarse_flags);

    # Solve coarse
    uH,UH=solve_maxwell(dmodel,order)
    check_error_maxwell(dmodel,order,uH)

    # Solve fine
    uh,Uh=solve_maxwell(fmodel,order)
    check_error_maxwell(fmodel,order,uh)

    # # prolongation via interpolation
    uHh=interpolate(uH,Uh)
    e = uh - uHh

    trian = Triangulation(fmodel)
    degree = 2*(order+1)
    dΩh = Measure(trian,degree)

    el2 = sqrt(sum( ∫( e⋅e )*dΩh ))
    tol=1e-6
    @assert el2 < tol
  end

  function test_2d(ranks,order)
    coarse_model=CartesianDiscreteModel((0,1,0,1),(1,1))
    dmodel=OctreeDistributedDiscreteModel(ranks,coarse_model,1)
    test_refine_and_coarsen_at_once(ranks,dmodel,order)
    rdmodel=dmodel
    for i=1:5
     rdmodel=test_transfer_ops_and_redistribute(ranks,rdmodel,order)
    end
  end 

  function test_3d(ranks,order)
    coarse_model=CartesianDiscreteModel((0,1,0,1,0,1),(2,2,2))
    dmodel=OctreeDistributedDiscreteModel(ranks,coarse_model,0)
    test_refine_and_coarsen_at_once(ranks,dmodel,order)
    rdmodel=dmodel
    for i=1:5
      rdmodel=test_transfer_ops_and_redistribute(ranks,rdmodel,order)
    end
  end

  u_ex_2D((x,y)) = 2*VectorValue(-y,x)
  f_ex_2D(x) = u_ex_2D(x)
  u_ex_3D((x,y,z)) = 2*VectorValue(-y,x,0.) - VectorValue(0.,-z,y)
  f_ex_3D(x) = u_ex_3D(x)
  
  function get_analytical_functions(Dc)
    if Dc==2
      return u_ex_2D, f_ex_2D
    else
      @assert Dc==3
      return u_ex_3D, f_ex_3D
    end
  end

  include("CoarseDiscreteModelsTools.jl")

  function solve_maxwell(model::GridapDistributed.DistributedDiscreteModel{Dc},order) where {Dc}
    u_ex, f_ex=get_analytical_functions(Dc)

    V = FESpace(model,
                ReferenceFE(nedelec,order),
                conformity=:Hcurl,
                dirichlet_tags="boundary")
    
    U = TrialFESpace(V,u_ex)
    
    trian = Triangulation(model)
    degree = 2*(order+1)
    dΩ = Measure(trian,degree)
        
    a(u,v) = ∫( (∇×u)⋅(∇×v) + u⋅v )dΩ
    l(v) = ∫(f_ex⋅v)dΩ

    op = AffineFEOperator(a,l,U,V)
    if (num_free_dofs(U)==0)
      # UMFPACK cannot handle empty linear systems
      uh = zero(U)
    else
      uh = solve(op)
    end
   
    # uh_ex=interpolate(u_ex_3D,U)
    # map(local_views(get_free_dof_values(uh_ex)), local_views(op.op.matrix), local_views(op.op.vector)) do U_ex, A, b 
    #   r_ex = A*U_ex - b
    #   println(r_ex)
    # end
    uh,U
  end 

  function check_error_maxwell(model::GridapDistributed.DistributedDiscreteModel{Dc},order,uh) where {Dc}
    trian = Triangulation(model)
    degree = 2*(order+1)
    dΩ = Measure(trian,degree)

    u_ex, f_ex = get_analytical_functions(Dc)
    
    eu = u_ex - uh

    l2(v) = sqrt(sum(∫(v⋅v)*dΩ))
    hcurl(v) = sqrt(sum(∫(v⋅v + (∇×v)⋅(∇×v))*dΩ))
    
    eu_l2 = l2(eu)
    eu_hcurl = hcurl(eu)
    
    tol = 1.0e-6
    @test eu_l2 < tol
    @test eu_hcurl < tol
  end 

  function run(distribute)
    # debug_logger = ConsoleLogger(stderr, Logging.Debug)
    # global_logger(debug_logger); # Enable the debug logger globally
    np    = MPI.Comm_size(MPI.COMM_WORLD)
    ranks = distribute(LinearIndices((np,)))

    for order=0:2
    test_2d(ranks,order)
    end

    for order=0:2
      test_3d(ranks,order)
    end
  end
end
