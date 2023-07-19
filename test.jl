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
ranks=distribute_with_mpi(LinearIndices((1,)))
coarse_model=CartesianDiscreteModel((0,1,0,1),(1,1))
dmodel=OctreeDistributedDiscreteModel(ranks,coarse_model,1)

map(ranks,partition(dmodel.dmodel.face_gids[end])) do rank, indices
    print("$(rank): $(local_to_owner(indices))"); print("\n") 
end 
reffe=ReferenceFE(lagrangian,Float64,order)
VH=FESpace(dmodel,reffe,conformity=:H1;dirichlet_tags="boundary")
UH=TrialFESpace(VH,u)
ref_coarse_flags=map(ranks,partition(get_cell_gids(dmodel.dmodel))) do rank,indices
    flags=zeros(Cint,length(indices))
    flags.=-1
    if (rank==1)
        #flags[1:2].=[refine_flag,nothing_flag]
        flags[1:4].=[refine_flag,nothing_flag,nothing_flag,nothing_flag]
    elseif (rank==2)
        flags[1:2].=[nothing_flag,nothing_flag]
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

rdmodel_red, red_glue=GridapDistributed.redistribute(fmodel);


# Vh=FESpace(fmodel,reffe,conformity=:H1;dirichlet_tags="boundary")
# map(ranks,partition(Vh.gids)) do rank, indices 
#     print("$(rank): $(local_to_owner(indices))"); print("\n")
#     print("$(rank): $(local_to_global(indices))"); print("\n")
# end 
# Uh=TrialFESpace(Vh,u)

# ΩH  = Triangulation(dmodel)
# dΩH = Measure(ΩH,degree)

# aH(u,v) = ∫( ∇(v)⊙∇(u) )*dΩH
# bH(v) = ∫(v*f)*dΩH

# op = AffineFEOperator(aH,bH,UH,VH)
# uH = solve(op)
# e = u - uH

# # # Compute errors
# el2 = sqrt(sum( ∫( e*e )*dΩH ))
# eh1 = sqrt(sum( ∫( e*e + ∇(e)⋅∇(e) )*dΩH ))

# tol=1e-8
# @assert el2 < tol
# @assert eh1 < tol


# Ωh  = Triangulation(fmodel)
# dΩh = Measure(Ωh,degree)

# ah(u,v) = ∫( ∇(v)⊙∇(u) )*dΩh
# bh(v) = ∫(v*f)*dΩh

# op = AffineFEOperator(ah,bh,Uh,Vh)
# uh = solve(op)
# e = u - uh

# # # Compute errors
# el2 = sqrt(sum( ∫( e*e )*dΩh ))
# eh1 = sqrt(sum( ∫( e*e + ∇(e)⋅∇(e) )*dΩh ))

# tol=1e-8
# @assert el2 < tol
# @assert eh1 < tol

# # prolongation via interpolation
# uHh=interpolate(uH,Uh)   
# e = uh - uHh
# el2 = sqrt(sum( ∫( e*e )*dΩh ))
# tol=1e-8
# @assert el2 < tol

# # prolongation via L2-projection 
# # Coarse FEFunction -> Fine FEFunction, by projection
# ah(u,v)  = ∫(v⋅u)*dΩh
# lh(v)    = ∫(v⋅uH)*dΩh
# oph      = AffineFEOperator(ah,lh,Uh,Vh)
# uHh      = solve(oph)
# e = uh - uHh
# el2 = sqrt(sum( ∫( e*e )*dΩh ))
# tol=1e-8
# @assert el2 < tol

# # restriction via interpolation
# uhH=interpolate(uh,UH) 
# e = uH - uhH
# el2 = sqrt(sum( ∫( e*e )*dΩh ))
# tol=1e-8
# @assert el2 < tol

# # restriction via L2-projection
# dΩhH = Measure(ΩH,Ωh,2*order)
# aH(u,v) = ∫(v⋅u)*dΩH
# lH(v)   = ∫(v⋅uh)*dΩhH
# oph     = AffineFEOperator(aH,lH,UH,VH)
# uhH     = solve(oph)
# e       = uH - uhH
# el2     = sqrt(sum( ∫( e*e )*dΩH ))



