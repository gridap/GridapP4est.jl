module NonConformingOctreeDistributedDiscreteModelsTests
  using P4est_wrapper
  using GridapP4est
  using Gridap
  using PartitionedArrays
  using GridapDistributed
  using MPI
  using Gridap.FESpaces
  using FillArrays
  using Logging

  function setup_model(::Type{Val{3}}, perm)
  #               5 +--------+ 7 
  #                /        /|
  #               /        / |
  #            6 +--------+  |
  #              |        |  |
  #              |  1     |  + 3 
  #              |        | /
  #              |        |/
  #            2 +--------+ 4

  #     6  +--------+ 8 
  #       /        /|
  #      /        / |
  #  11 +--------+  |
  #     |        |  |
  #     |  2     |  + 4
  #     |        | /
  #     |        |/
  #   9 +--------+ 10
    ptr  = [ 1, 9, 17 ]
    if (perm==1)
      data = [ 1,2,3,4,5,6,7,8, 2,9,4,10,6,11,8,12 ]
    elseif (perm==2)
      data = [ 1,2,3,4,5,6,7,8, 10,12,4,8,9,11,2,6 ]
    elseif (perm==3)
      data = [ 1,2,3,4,5,6,7,8, 12,11,8,6,10,9,4,2 ]
    elseif (perm==4) 
      data = [ 1,2,3,4,5,6,7,8, 11,9,6,2,12,10,8,4 ]
    end  
    cell_vertex_lids = Gridap.Arrays.Table(data,ptr)
    node_coordinates = Vector{Point{3,Float64}}(undef,12)
    node_coordinates[1]=Point{3,Float64}(0.0,0.0,0.0)
    node_coordinates[2]=Point{3,Float64}(1.0,0.0,0.0)
    node_coordinates[3]=Point{3,Float64}(0.0,1.0,0.0)
    node_coordinates[4]=Point{3,Float64}(1.0,1.0,0.0)
    node_coordinates[5]=Point{3,Float64}(0.0,0.0,1.0)
    node_coordinates[6]=Point{3,Float64}(1.0,0.0,1.0)
    node_coordinates[7]=Point{3,Float64}(0.0,1.0,1.0)
    node_coordinates[8]=Point{3,Float64}(1.0,1.0,1.0)
    node_coordinates[9]=Point{3,Float64}(2.0,0.0,0.0)
    node_coordinates[10]=Point{3,Float64}(2.0,1.0,0.0)
    node_coordinates[11]=Point{3,Float64}(2.0,0.0,1.0)
    node_coordinates[12]=Point{3,Float64}(2.0,1.0,1.0)

    polytope=HEX
    scalar_reffe=Gridap.ReferenceFEs.ReferenceFE(polytope,Gridap.ReferenceFEs.lagrangian,Float64,1)
    cell_types=collect(Fill(1,length(cell_vertex_lids)))
    cell_reffes=[scalar_reffe]
    grid = Gridap.Geometry.UnstructuredGrid(node_coordinates,
                                            cell_vertex_lids,
                                            cell_reffes,
                                            cell_types,
                                            Gridap.Geometry.NonOriented())
    m=Gridap.Geometry.UnstructuredDiscreteModel(grid)
    labels = get_face_labeling(m)
    labels.d_to_dface_to_entity[1].=2
    if (perm==1 || perm==2)
      labels.d_to_dface_to_entity[2].=2
      labels.d_to_dface_to_entity[3].=[2,2,2,2,2,1,2,2,2,2,2]
    elseif (perm==3 || perm==4)
      labels.d_to_dface_to_entity[2].=2
      labels.d_to_dface_to_entity[3].=[2,2,2,2,2,1,2,2,2,2,2]
    end
    labels.d_to_dface_to_entity[4].=1  
    add_tag!(labels,"boundary",[2])
    add_tag!(labels,"interior",[1])
    m
  end

  function setup_model(::Type{Val{2}}, perm)
    @assert perm ∈ (1,2,3,4)
    #
    #  3-------4-------6
    #  |       |       |
    #  |       |       |
    #  |       |       |
    #  1-------2-------5
    #
        ptr  = [ 1, 5, 9 ]
        if (perm==1)
          data = [ 1,2,3,4, 2,5,4,6 ]
        elseif (perm==2)
          data = [ 1,2,3,4, 6,4,5,2 ]
        elseif (perm==3)
          data = [ 4,3,2,1, 2,5,4,6 ]
        elseif (perm==4) 
          data = [ 4,3,2,1, 6,4,5,2 ]
        end  
        cell_vertex_lids = Gridap.Arrays.Table(data,ptr)
        node_coordinates = Vector{Point{2,Float64}}(undef,6)
        node_coordinates[1]=Point{2,Float64}(0.0,0.0)
        node_coordinates[2]=Point{2,Float64}(1.0,0.0)
        node_coordinates[3]=Point{2,Float64}(0.0,1.0)
        node_coordinates[4]=Point{2,Float64}(1.0,1.0)
        node_coordinates[5]=Point{2,Float64}(2.0,0.0)
        node_coordinates[6]=Point{2,Float64}(2.0,1.0)
    
        polytope=QUAD
        scalar_reffe=Gridap.ReferenceFEs.ReferenceFE(polytope,Gridap.ReferenceFEs.lagrangian,Float64,1)
        cell_types=collect(Fill(1,length(cell_vertex_lids)))
        cell_reffes=[scalar_reffe]
        grid = Gridap.Geometry.UnstructuredGrid(node_coordinates,
                                                cell_vertex_lids,
                                                cell_reffes,
                                                cell_types,
                                                Gridap.Geometry.NonOriented())
        m=Gridap.Geometry.UnstructuredDiscreteModel(grid)
        labels = get_face_labeling(m)
        labels.d_to_dface_to_entity[1].=2
        if (perm==1 || perm==2)
          labels.d_to_dface_to_entity[2].=[2,2,2,1,2,2,2]
        elseif (perm==3 || perm==4)
          labels.d_to_dface_to_entity[2].=[2,2,1,2,2,2,2] 
        end
        labels.d_to_dface_to_entity[3].=1
        add_tag!(labels,"boundary",[2])
        add_tag!(labels,"interior",[1])
        m
  end


  function test_transfer_ops_and_redistribute(ranks,dmodel,order)
    # Define manufactured functions
    u(x) = x[1]+x[2]^order
    f(x) = -Δ(u)(x)
    degree = 2*order+1
    reffe=ReferenceFE(lagrangian,Float64,order)
    VH=FESpace(dmodel,reffe,conformity=:H1;dirichlet_tags="boundary")
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
    fmodel,glue=adapt(dmodel,ref_coarse_flags);
    # map(ranks,glue) do rank, glue 
    #   if rank==2
    #     print(glue.n2o_faces_map[end]); print("\n")
    #   end  
    # end
    Vh=FESpace(fmodel,reffe,conformity=:H1;dirichlet_tags="boundary")
    Uh=TrialFESpace(Vh,u)

    ΩH  = Triangulation(dmodel)
    dΩH = Measure(ΩH,degree)

    aH(u,v) = ∫( ∇(v)⊙∇(u) )*dΩH
    bH(v) = ∫(v*f)*dΩH

    op = AffineFEOperator(aH,bH,UH,VH)
    uH = solve(op)
    e = u - uH

    # # Compute errors
    el2 = sqrt(sum( ∫( e*e )*dΩH ))
    eh1 = sqrt(sum( ∫( e*e + ∇(e)⋅∇(e) )*dΩH ))

    tol=1e-8
    @assert el2 < tol
    @assert eh1 < tol


    Ωh  = Triangulation(fmodel)
    dΩh = Measure(Ωh,degree)

    ah(u,v) = ∫( ∇(v)⊙∇(u) )*dΩh
    bh(v) = ∫(v*f)*dΩh

    op = AffineFEOperator(ah,bh,Uh,Vh)
    uh = solve(op)
    e = u - uh

    # # Compute errors
    el2 = sqrt(sum( ∫( e*e )*dΩh ))
    eh1 = sqrt(sum( ∫( e*e + ∇(e)⋅∇(e) )*dΩh ))

    tol=1e-8
    @assert el2 < tol
    @assert eh1 < tol

    # prolongation via interpolation
    uHh=interpolate(uH,Uh)   
    e = uh - uHh
    el2 = sqrt(sum( ∫( e*e )*dΩh ))
    tol=1e-8
    @assert el2 < tol

    # prolongation via L2-projection 
    # Coarse FEFunction -> Fine FEFunction, by projection
    ahp(u,v)  = ∫(v⋅u)*dΩh
    lhp(v)    = ∫(v⋅uH)*dΩh
    oph      = AffineFEOperator(ahp,lhp,Uh,Vh)
    uHh      = solve(oph)
    e = uh - uHh
    el2 = sqrt(sum( ∫( e*e )*dΩh ))
    tol=1e-8
    @assert el2 < tol

    # restriction via interpolation
    uhH=interpolate(uh,UH) 
    e = uH - uhH
    el2 = sqrt(sum( ∫( e*e )*dΩh ))
    tol=1e-8
    @assert el2 < tol

    # restriction via L2-projection
    dΩhH = Measure(ΩH,Ωh,2*order)
    aHp(u,v) = ∫(v⋅u)*dΩH
    lHp(v)   = ∫(v⋅uh)*dΩhH
    oph     = AffineFEOperator(aHp,lHp,UH,VH)
    uhH     = solve(oph)
    e       = uH - uhH
    el2     = sqrt(sum( ∫( e*e )*dΩH ))

    fmodel_red, red_glue=GridapDistributed.redistribute(fmodel);
    Vhred=FESpace(fmodel_red,reffe,conformity=:H1;dirichlet_tags="boundary")
    Uhred=TrialFESpace(Vhred,u)

    Ωhred  = Triangulation(fmodel_red)
    dΩhred = Measure(Ωhred,degree)

    ahred(u,v) = ∫( ∇(v)⊙∇(u) )*dΩhred
    bhred(v)   = ∫(v*f)*dΩhred

    op    = AffineFEOperator(ahred,bhred,Uhred,Vhred)
    uhred = solve(op)
    e = u - uhred
    el2 = sqrt(sum( ∫( e*e )*dΩhred ))
    @assert el2 < tol


    uhred2 = GridapP4est.redistribute_fe_function(uh,Vhred,fmodel_red,red_glue)
    e = u - uhred2
    el2 = sqrt(sum( ∫( e*e )*dΩhred ))
    tol=1e-8
    @assert el2 < tol

    fmodel_red
  end

  function test_refine_and_coarsen_at_once(ranks,
            dmodel::OctreeDistributedDiscreteModel{Dc},
            order) where Dc
    u(x) = x[1]+x[2]^order
    f(x) = -Δ(u)(x)
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
    fmodel,glue=adapt(dmodel,ref_coarse_flags);

    reffe=ReferenceFE(lagrangian,Float64,order)
    VH=FESpace(dmodel,reffe,conformity=:H1;dirichlet_tags="boundary")
    UH=TrialFESpace(VH,u)

    Vh=FESpace(fmodel,reffe,conformity=:H1;dirichlet_tags="boundary")
    Uh=TrialFESpace(Vh,u)
    ΩH  = Triangulation(dmodel)
    dΩH = Measure(ΩH,degree)

    aH(u,v) = ∫( ∇(v)⊙∇(u) )*dΩH
    bH(v) = ∫(v*f)*dΩH

    op = AffineFEOperator(aH,bH,UH,VH)
    uH = solve(op)
    e = u - uH

    # # Compute errors
    el2 = sqrt(sum( ∫( e*e )*dΩH ))
    eh1 = sqrt(sum( ∫( e*e + ∇(e)⋅∇(e) )*dΩH ))

    tol=1e-8
    @assert el2 < tol
    @assert eh1 < tol

    Ωh  = Triangulation(fmodel)
    dΩh = Measure(Ωh,degree)

    ah(u,v) = ∫( ∇(v)⊙∇(u) )*dΩh
    bh(v) = ∫(v*f)*dΩh

    op = AffineFEOperator(ah,bh,Uh,Vh)
    uh = solve(op)
    e = u - uh

    # # Compute errors
    el2 = sqrt(sum( ∫( e*e )*dΩh ))
    eh1 = sqrt(sum( ∫( e*e + ∇(e)⋅∇(e) )*dΩh ))

    tol=1e-8
    @assert el2 < tol
    @assert eh1 < tol

    # prolongation via interpolation
    uHh=interpolate(uH,Uh)
    e = uh - uHh
    el2 = sqrt(sum( ∫( e*e )*dΩh ))
    tol=1e-8
    @assert el2 < tol
  end

  function test_2d(ranks,order)
    coarse_model=CartesianDiscreteModel((0,1,0,1),(1,1))
    dmodel=OctreeDistributedDiscreteModel(ranks,coarse_model,2)
    test_refine_and_coarsen_at_once(ranks,dmodel,order)
    rdmodel=dmodel
    # for i=1:5
    #  rdmodel=test_transfer_ops_and_redistribute(ranks,rdmodel,order)
    # end
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

  function test(ranks,TVDc::Type{Val{Dc}}, perm, order) where Dc
    # function test_solve(dmodel,order)
    #   # Define manufactured functions
    #   u(x) = x[1]+x[2]^order
    #   f(x) = -Δ(u)(x)

    #   # FE Spaces
    #   reffe = ReferenceFE(lagrangian,Float64,order)
    #   V = TestFESpace(dmodel,reffe,dirichlet_tags="boundary")
    #   U = TrialFESpace(V,u)

    #   # Define integration mesh and quadrature
    #   degree = 2*order+1
    #   Ω = Triangulation(dmodel)
    #   dΩ = Measure(Ω,degree)

    #   a(u,v) = ∫( ∇(v)⊙∇(u) )*dΩ
    #   b(v) = ∫(v*f)*dΩ

    #   op = AffineFEOperator(a,b,U,V)
    #   uh = solve(op)

    #   e = u - uh

    #   # Compute errors
    #   el2 = sqrt(sum( ∫( e*e )*dΩ ))
    #   eh1 = sqrt(sum( ∫( e*e + ∇(e)⋅∇(e) )*dΩ ))

    #   tol=1e-6
    #   println("$(el2) < $(tol)")
    #   println("$(eh1) < $(tol)")
    #   @assert el2 < tol
    #   @assert eh1 < tol
    # end

    coarse_model = setup_model(TVDc,perm)
    model = OctreeDistributedDiscreteModel(ranks, coarse_model, 0)
    # test_transfer_ops_and_redistribute(ranks,model,order)




    #test_solve(model,order)
    # ref_coarse_flags=map(ranks) do _
    #   [refine_flag,nothing_flag]
    # end 
    # dmodel,adaptivity_glue=adapt(model,ref_coarse_flags)
    # non_conforming_glue=dmodel.non_conforming_glue
    # test_solve(dmodel,order)
    # for i=1:3
    #   ref_coarse_flags=map(ranks,partition(get_cell_gids(dmodel.dmodel))) do rank,indices
    #     flags=zeros(Cint,length(indices))
    #     flags.=nothing_flag
    #     flags[1]=refine_flag
    #     flags[own_length(indices)]=refine_flag    
    #     print("rank: $(rank) flags: $(flags)"); print("\n")
    #     flags
    #   end 
    #   dmodel,glue=adapt(dmodel,ref_coarse_flags);
    #   # test_solve(dmodel,order)
    # end



  end

  function run(distribute)
    debug_logger = ConsoleLogger(stderr, Logging.Debug)
    global_logger(debug_logger); # Enable the debug logger globally

    ranks = distribute(LinearIndices((MPI.Comm_size(MPI.COMM_WORLD),)))
    # for Dc=2:3, perm=1:4, order=1:4
    #    test(ranks,Val{Dc},perm,order)
    # end
    for order=1:1
      test_2d(ranks,order)
      #test_3d(ranks,order)
    end


  end

end