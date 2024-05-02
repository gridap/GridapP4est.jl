module DarcyNonConformingOctreeModelsTests
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
    xH,XH=solve_darcy(dmodel,order)
    check_error_darcy(dmodel,order,xH)
    uH,pH=xH
    UH,PH=XH

    # Solve fine
    xh,Xh=solve_darcy(fmodel,order)
    check_error_darcy(fmodel,order,xh)
    uh,ph=xh
    Uh,Ph=Xh

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
    xh_red,Xh_red=solve_darcy(fmodel_red,order)
    check_error_darcy(fmodel_red,order,xh_red)
    uhred,phred=xh_red
    Uh_red,Ph_red=Xh_red


    trian = Triangulation(fmodel_red)
    degree = 2*(order+1)
    dΩhred = Measure(trian,degree)

    u_ex, p_ex, f_ex=get_analytical_functions(Dc)

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
           flags[1:2^Dc].=coarsen_flag
        else 
           flags[own_length(indices)]=refine_flag
        end 
        flags
    end
    fmodel,glue=Gridap.Adaptivity.adapt(dmodel,ref_coarse_flags);

    # Solve coarse
    xH,XH=solve_darcy(dmodel,order)
    check_error_darcy(dmodel,order,xH)
    uH,pH=xH

    # Solve fine
    xh,Xh=solve_darcy(fmodel,order)
    check_error_darcy(fmodel,order,xh)
    uh,ph=xh

    # prolongation via interpolation
    Uh,Ph=Xh
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
    dmodel=OctreeDistributedDiscreteModel(ranks,coarse_model,2)
    test_refine_and_coarsen_at_once(ranks,dmodel,order)
    rdmodel=dmodel
    for i=1:5
     rdmodel=test_transfer_ops_and_redistribute(ranks,rdmodel,order)
    end
  end 

  function test_3d(ranks,order)
    coarse_model=CartesianDiscreteModel((0,1,0,1,0,1),(1,1,1))
    dmodel=OctreeDistributedDiscreteModel(ranks,coarse_model,2)
    test_refine_and_coarsen_at_once(ranks,dmodel,order)
    rdmodel=dmodel
    for i=1:5
      rdmodel=test_transfer_ops_and_redistribute(ranks,rdmodel,order)
    end
  end

  u_ex_2D(x) = VectorValue(2*x[1],x[1]+x[2])
  p_ex_2D(x) = x[1]-x[2]
  f_ex_2D(x) = u_ex_2D(x) + ∇(p_ex_2D)(x)
  u_ex_3D(x) = VectorValue(x[1],x[2],x[3])
  p_ex_3D(x) = x[2]
  f_ex_3D(x) = u_ex_3D(x) + ∇(p_ex_3D)(x)
  
  function get_analytical_functions(Dc)
    if Dc==2
      return u_ex_2D, p_ex_2D, f_ex_2D
    else
      @assert Dc==3
      return u_ex_3D, p_ex_3D, f_ex_3D
    end
  end

  include("CoarseDiscreteModelsTools.jl")

  function solve_darcy(model::GridapDistributed.DistributedDiscreteModel{Dc},order) where {Dc}
    if (Dc==2)
      dirichlet_tags=[5,6]
      neumanntags = [7,8]
    else
      @assert Dc==3
      dirichlet_tags=[21,22,23]
      neumanntags=[24,25,26]
    end 

    u_ex, p_ex, f_ex=get_analytical_functions(Dc)

    V = FESpace(model,
                ReferenceFE(raviart_thomas,Float64,order),
                conformity=:Hdiv,
                dirichlet_tags=dirichlet_tags)
    
    Q = FESpace(model,
                ReferenceFE(lagrangian,Float64,order); 
                conformity=:L2)
    
    U = TrialFESpace(V,u_ex)
    P = TrialFESpace(Q)
    
    Y = MultiFieldFESpace([V, Q])
    X = MultiFieldFESpace([U, P])
    
    trian = Triangulation(model)
    degree = 2*(order+1)
    dΩ = Measure(trian,degree)
    
    btrian = Boundary(model,tags=neumanntags)
    degree = 2*(order+1)
    dΓ = Measure(btrian,degree)
    nb = get_normal_vector(btrian)
    
    a((u, p),(v, q)) = ∫(u⋅v)dΩ +∫(q*(∇⋅u))dΩ-∫((∇⋅v)*p)dΩ
    b(( v, q)) = ∫( v⋅f_ex + q*(∇⋅u_ex))dΩ - ∫((v⋅nb)*p_ex )dΓ

    op = AffineFEOperator(a,b,X,Y)
    xh = solve(op)
    xh,X
  end 

  function check_error_darcy(model::GridapDistributed.DistributedDiscreteModel{Dc},order,xh) where {Dc}
    trian = Triangulation(model)
    degree = 2*(order+1)
    dΩ = Measure(trian,degree)

    uh, ph = xh

    u_ex, p_ex, f_ex = get_analytical_functions(Dc)
    
    eu = u_ex - uh
    ep = p_ex - ph

    l2(v) = sqrt(sum(∫(v⋅v)*dΩ))
    h1(v) = sqrt(sum(∫(v*v + ∇(v)⋅∇(v))*dΩ))
    
    eu_l2 = l2(eu)
    ep_l2 = l2(ep)
    ep_h1 = h1(ep)
    
    tol = 1.0e-6
    @test eu_l2 < tol
    @test ep_l2 < tol
    @test ep_h1 < tol
  end

  function run(distribute)
    # debug_logger = ConsoleLogger(stderr, Logging.Debug)
    # global_logger(debug_logger); # Enable the debug logger globally
    np    = MPI.Comm_size(MPI.COMM_WORLD)
    ranks = distribute(LinearIndices((np,)))

    order=1
    test_2d(ranks,order)

    order=2
    test_3d(ranks,order)
  end
end