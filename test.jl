using PartitionedArrays
using Gridap
using GridapP4est 
using GridapDistributed
using MPI


# Define integration mesh and quadrature
order=1
# Define manufactured functions
u(x) = x[1]+x[2]^order
f(x) = -Δ(u)(x)
degree = 2*order+1

MPI.Init()
ranks=distribute_with_mpi(LinearIndices((2,)))
coarse_model=CartesianDiscreteModel((0,1,0,1),(1,1))
dmodel=OctreeDistributedDiscreteModel(ranks,coarse_model,2)

function test(ranks,dmodel)
    map(ranks,partition(dmodel.dmodel.face_gids[end])) do rank, indices
        print("$(rank): $(local_to_owner(indices))"); print("\n") 
    end 
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

        print("rank: $(rank) flags: $(flags)"); print("\n")
        flags
    end 
    fmodel,glue=refine(dmodel,ref_coarse_flags);
    map(ranks,glue) do rank, glue 
        if rank==2
        print(glue.n2o_faces_map[end]); print("\n")
        print(glue.refinement_rules); print("\n")
        end  
    end

    Vh=FESpace(fmodel,reffe,conformity=:H1;dirichlet_tags="boundary")
    map(ranks,partition(Vh.gids)) do rank, indices 
        print("$(rank): $(local_to_owner(indices))"); print("\n")
        print("$(rank): $(local_to_global(indices))"); print("\n")
    end 
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


    map(Vhred.spaces) do space
        if (MPI.Comm_rank(MPI.COMM_WORLD)==1)
           println("XXXX: $(get_cell_dof_ids(space.space))")
        end  
    end 

    uhred2 = GridapP4est.redistribute_fe_function(uh,Vhred,fmodel_red,red_glue)
    e = u - uhred2
    el2 = sqrt(sum( ∫( e*e )*dΩhred ))
    tol=1e-8
    println(el2)
    @assert el2 < tol

    fmodel_red
end

function test_coarsen(ranks,dmodel)
    ref_coarse_flags=map(ranks,partition(get_cell_gids(dmodel.dmodel))) do rank,indices
        flags=zeros(Cint,length(indices))
        flags.=nothing_flag        
        if (rank==1)
           flags[1:4].=coarsen_flag
        else 
           flags[own_length(indices)]=refine_flag
        end 
        flags
    end
    fmodel,glue=refine(dmodel,ref_coarse_flags);

    order=1
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
test_coarsen(ranks,dmodel);

function f()
  rdmodel=dmodel
  for i=1:5
   rdmodel=test(ranks,rdmodel)
  end
end
f()


