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

#   function test_transfer_ops_and_redistribute(ranks,dmodel,order)
#     # Define manufactured functions
#     u(x) = x[1]+x[2]^order
#     f(x) = -Δ(u)(x)
#     degree = 2*order+1
#     reffe=ReferenceFE(lagrangian,Float64,order)
#     VH=FESpace(dmodel,reffe,conformity=:H1;dirichlet_tags="boundary")
#     UH=TrialFESpace(VH,u)
#     ref_coarse_flags=map(ranks,partition(get_cell_gids(dmodel.dmodel))) do rank,indices
#         flags=zeros(Cint,length(indices))
#         flags.=nothing_flag
        
#         flags[1]=refine_flag
#         flags[own_length(indices)]=refine_flag

#         # To create some unbalance
#         if (rank%2==0 && own_length(indices)>1)
#             flags[div(own_length(indices),2)]=refine_flag
#         end
#         flags
#     end 
#     fmodel,glue=adapt(dmodel,ref_coarse_flags);
#     # map(ranks,glue) do rank, glue 
#     #   if rank==2
#     #     print(glue.n2o_faces_map[end]); print("\n")
#     #   end  
#     # end
#     Vh=FESpace(fmodel,reffe,conformity=:H1;dirichlet_tags="boundary")
#     Uh=TrialFESpace(Vh,u)

#     ΩH  = Triangulation(dmodel)
#     dΩH = Measure(ΩH,degree)

#     aH(u,v) = ∫( ∇(v)⊙∇(u) )*dΩH
#     bH(v) = ∫(v*f)*dΩH

#     op = AffineFEOperator(aH,bH,UH,VH)
#     uH = solve(op)
#     e = u - uH

#     # # Compute errors
#     el2 = sqrt(sum( ∫( e*e )*dΩH ))
#     eh1 = sqrt(sum( ∫( e*e + ∇(e)⋅∇(e) )*dΩH ))

#     tol=1e-6
#     println("[SOLVE COARSE] el2 < tol: $(el2) < $(tol)")
#     println("[SOLVE COARSE] eh1 < tol: $(eh1) < $(tol)")
#     @assert el2 < tol
#     @assert eh1 < tol


#     Ωh  = Triangulation(fmodel)
#     dΩh = Measure(Ωh,degree)

#     ah(u,v) = ∫( ∇(v)⊙∇(u) )*dΩh
#     bh(v) = ∫(v*f)*dΩh

#     op = AffineFEOperator(ah,bh,Uh,Vh)
#     uh = solve(op)
#     e = u - uh

#     writevtk(ΩH, "ctrian", cellfields=["uH"=>uH])
#     writevtk(Ωh, "ftrian", cellfields=["uh"=>uh])

#     # # Compute errors

#     el2 = sqrt(sum( ∫( e*e )*dΩh ))
#     eh1 = sqrt(sum( ∫( e*e + ∇(e)⋅∇(e) )*dΩh ))
 
#     println("[SOLVE FINE] el2 < tol: $(el2) < $(tol)")
#     println("[SOLVE FINE] eh1 < tol: $(eh1) < $(tol)")
#     @assert el2 < tol
#     @assert eh1 < tol

#     # prolongation via interpolation
#     uHh=interpolate(uH,Uh)   
#     e = uh - uHh
#     el2 = sqrt(sum( ∫( e*e )*dΩh ))
#     tol=1e-6
#     println("[INTERPOLATION] el2 < tol: $(el2) < $(tol)")
#     @assert el2 < tol

#     # prolongation via L2-projection 
#     # Coarse FEFunction -> Fine FEFunction, by projection
#     ahp(u,v)  = ∫(v⋅u)*dΩh
#     lhp(v)    = ∫(v⋅uH)*dΩh
#     oph      = AffineFEOperator(ahp,lhp,Uh,Vh)
#     uHh      = solve(oph)
#     e = uh - uHh
#     el2 = sqrt(sum( ∫( e*e )*dΩh ))
#     println("[L2 PROJECTION] el2 < tol: $(el2) < $(tol)")
#     @assert el2 < tol

#     # restriction via interpolation
#     uhH=interpolate(uh,UH) 
#     e = uH - uhH
#     el2 = sqrt(sum( ∫( e*e )*dΩh ))
#     println("[INTERPOLATION] el2 < tol: $(el2) < $(tol)")
#     @assert el2 < tol

#     # restriction via L2-projection
#     dΩhH = Measure(ΩH,Ωh,2*order)
#     aHp(u,v) = ∫(v⋅u)*dΩH
#     lHp(v)   = ∫(v⋅uh)*dΩhH
#     oph     = AffineFEOperator(aHp,lHp,UH,VH)
#     uhH     = solve(oph)
#     e       = uH - uhH
#     el2     = sqrt(sum( ∫( e*e )*dΩH ))

#     fmodel_red, red_glue=GridapDistributed.redistribute(fmodel);
#     Vhred=FESpace(fmodel_red,reffe,conformity=:H1;dirichlet_tags="boundary")
#     Uhred=TrialFESpace(Vhred,u)

#     Ωhred  = Triangulation(fmodel_red)
#     dΩhred = Measure(Ωhred,degree)

#     ahred(u,v) = ∫( ∇(v)⊙∇(u) )*dΩhred
#     bhred(v)   = ∫(v*f)*dΩhred

#     op    = AffineFEOperator(ahred,bhred,Uhred,Vhred)
#     uhred = solve(op)
#     e = u - uhred
#     el2 = sqrt(sum( ∫( e*e )*dΩhred ))
#     println("[SOLVE FINE REDISTRIBUTED] el2 < tol: $(el2) < $(tol)")
#     @assert el2 < tol


#     uhred2 = GridapDistributed.redistribute_fe_function(uh,Vhred,fmodel_red,red_glue)
#     e = u - uhred2
#     el2 = sqrt(sum( ∫( e*e )*dΩhred ))
#     println("[REDISTRIBUTE SOLUTION] el2 < tol: $(el2) < $(tol)")
#     @assert el2 < tol

#     fmodel_red
#   end

#   function test_refine_and_coarsen_at_once(ranks,
#             dmodel::OctreeDistributedDiscreteModel{Dc},
#             order) where Dc
#     u(x) = x[1]+x[2]^order
#     f(x) = -Δ(u)(x)
#     degree = 2*order+1
#     ref_coarse_flags=map(ranks,partition(get_cell_gids(dmodel.dmodel))) do rank,indices
#         flags=zeros(Cint,length(indices))
#         flags.=nothing_flag        
#         if (rank==1)
#            flags[1:2^Dc].=coarsen_flag
#         else 
#            flags[own_length(indices)]=refine_flag
#         end 
#         flags
#     end
#     fmodel,glue=adapt(dmodel,ref_coarse_flags);

#     reffe=ReferenceFE(lagrangian,Float64,order)
#     VH=FESpace(dmodel,reffe,conformity=:H1;dirichlet_tags="boundary")
#     UH=TrialFESpace(VH,u)

#     Vh=FESpace(fmodel,reffe,conformity=:H1;dirichlet_tags="boundary")
#     Uh=TrialFESpace(Vh,u)
#     ΩH  = Triangulation(dmodel)
#     dΩH = Measure(ΩH,degree)

#     aH(u,v) = ∫( ∇(v)⊙∇(u) )*dΩH
#     bH(v) = ∫(v*f)*dΩH

#     op = AffineFEOperator(aH,bH,UH,VH)
#     uH = solve(op)
#     e = u - uH

#     # # Compute errors
#     el2 = sqrt(sum( ∫( e*e )*dΩH ))
#     eh1 = sqrt(sum( ∫( e*e + ∇(e)⋅∇(e) )*dΩH ))

#     tol=1e-7
#     @assert el2 < tol
#     @assert eh1 < tol

#     Ωh  = Triangulation(fmodel)
#     dΩh = Measure(Ωh,degree)

#     ah(u,v) = ∫( ∇(v)⊙∇(u) )*dΩh
#     bh(v) = ∫(v*f)*dΩh

#     op = AffineFEOperator(ah,bh,Uh,Vh)
#     uh = solve(op)
#     e = u - uh

#     # # Compute errors
#     el2 = sqrt(sum( ∫( e*e )*dΩh ))
#     eh1 = sqrt(sum( ∫( e*e + ∇(e)⋅∇(e) )*dΩh ))

#     tol=1e-7
#     @assert el2 < tol
#     @assert eh1 < tol

#     # prolongation via interpolation
#     uHh=interpolate(uH,Uh)
#     e = uh - uHh
#     el2 = sqrt(sum( ∫( e*e )*dΩh ))
#     tol=1e-7
#     @assert el2 < tol
#   end

  function test_2d(ranks,order)
    coarse_model=CartesianDiscreteModel((0,1,0,1),(1,1))
    dmodel=OctreeDistributedDiscreteModel(ranks,coarse_model,2)
    test_refine_and_coarsen_at_once(ranks,dmodel,order)
    rdmodel=dmodel
    for i=1:5
     rdmodel=test_transfer_ops_and_redistribute(ranks,rdmodel,order)
    end
  end 

#   function test_3d(ranks,order)
#     coarse_model=CartesianDiscreteModel((0,1,0,1,0,1),(1,1,1))
#     dmodel=OctreeDistributedDiscreteModel(ranks,coarse_model,2)
#     test_refine_and_coarsen_at_once(ranks,dmodel,order)
#     rdmodel=dmodel
#     for i=1:5
#       rdmodel=test_transfer_ops_and_redistribute(ranks,rdmodel,order)
#     end
#   end 

#   function test(ranks,TVDc::Type{Val{Dc}}, perm, order) where Dc
#     coarse_model = setup_model(TVDc,perm)
#     model = OctreeDistributedDiscreteModel(ranks, coarse_model, 1)
#     test_transfer_ops_and_redistribute(ranks,model,order)
#   end

  u_ex_2D(x) = VectorValue(2*x[1],x[1]+x[2])
  p_ex_2D(x) = x[1]-x[2]
  f_ex_2D(x) = u_ex_2D(x) + ∇(p_ex_2D)(x)
  u_ex_3D(x) = VectorValue(2*x[1],x[1]+x[2],4*x[3])
  p_ex_3D(x) = x[1]-x[2]+2x[3]
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

  function GridapDistributed.remove_ghost_cells(
    trian::Gridap.Adaptivity.AdaptedTriangulation{Dc,Dp,<:Union{SkeletonTriangulation,BoundaryTriangulation}},gids) where {Dc,Dp}
    GridapDistributed.remove_ghost_cells(trian.trian,gids)
  end

  function solve_darcy(model::GridapDistributed.DistributedDiscreteModel{Dc},order) where {Dc}
    if (Dc==2)
      dirichlet_tags=[5,6]
      neumanntags = [7,8]
    else
      @assert Dc==3
      dirichlet_tags=[21,22,23]
      neumanntags = [24,25,26]
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
  end 

  function check_error_darcy(model::GridapDistributed.DistributedDiscreteModel{Dc},order,xh) where {Dc}
    trian = Triangulation(model)
    degree = 2*(order+1)
    dΩ = Measure(trian,degree)

    uh = xh[1]
    uh, ph = xh

    u_ex, p_ex, f_ex = get_analytical_functions(Dc)
    
    eu = u_ex - uh
    ep = p_ex - ph

    l2(v) = sqrt(sum(∫(v⋅v)*dΩ))
    h1(v) = sqrt(sum(∫(v*v + ∇(v)⋅∇(v))*dΩ))
    
    eu_l2 = l2(eu)
    ep_l2 = l2(ep)
    ep_h1 = h1(ep)
    
    tol = 1.0e-9
    @test eu_l2 < tol
    @test ep_l2 < tol
    @test ep_h1 < tol
  end 

  function driver(ranks,coarse_model,order)
    model=OctreeDistributedDiscreteModel(ranks,coarse_model,1)
    ref_coarse_flags=map(ranks,partition(get_cell_gids(model.dmodel))) do rank,indices
        flags=zeros(Cint,length(indices))
        flags.=nothing_flag        
        #flags[1:2^2].=coarsen_flag
        flags[own_length(indices)]=refine_flag 
        flags
    end
    fmodel,glue=adapt(model,ref_coarse_flags)
    xh = solve_darcy(fmodel,order)
    check_error_darcy(fmodel,order,xh)
  end 

  function run(distribute)
    # debug_logger = ConsoleLogger(stderr, Logging.Debug)
    # global_logger(debug_logger); # Enable the debug logger globally
    np    = MPI.Comm_size(MPI.COMM_WORLD)
    ranks = distribute(LinearIndices((np,)))
    
    # domain = (0,1,0,1)
    # order = 1
    # coarse_model=CartesianDiscreteModel(domain,(1,1))
    # driver(ranks,coarse_model,order)

    domain = (0,1,0,1,0,1)
    order = 1 
    coarse_model=CartesianDiscreteModel(domain,(1,1,1))
    driver(ranks,coarse_model,order)
  end
end