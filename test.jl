using PartitionedArrays
using Gridap
using GridapP4est 
using GridapDistributed
using MPI

#MPI.Init()
ranks=distribute_with_mpi(LinearIndices((2,)))
coarse_model=CartesianDiscreteModel((0,1,0,1),(1,1))
dmodel=OctreeDistributedDiscreteModel(ranks,coarse_model,1)

map(ranks,partition(dmodel.dmodel.face_gids[end])) do rank, indices
    print("$(rank): $(local_to_owner(indices))"); print("\n") 
end 
reffe=ReferenceFE(lagrangian,Float64,1)
UH=FESpace(dmodel,reffe,conformity=:H1;dirichlet_tags=Int[])
ref_coarse_flags=map(ranks,partition(get_cell_gids(dmodel.dmodel))) do rank,indices
    flags=zeros(Cint,length(indices))
    flags.=-1
    if (rank==1)
        flags[1:2].=[refine_flag,nothing_flag]
    elseif (rank==2)
        flags[1:2].=[nothing_flag,refine_flag]
    end
    print("rank: $(rank) flags: $(flags)"); print("\n")
    flags
end 
rdmodel,glue=refine(dmodel,ref_coarse_flags);
map(ranks,glue) do rank, glue 
    if rank==2
      print(glue.n2o_faces_map[end]); print("\n")
      print(glue.refinement_rules); print("\n")
    end  
end
Vh=FESpace(rdmodel,reffe,conformity=:H1;dirichlet_tags=Int[])
map(ranks,partition(Vh.gids)) do rank, indices 
    print("$(rank): $(local_to_owner(indices))"); print("\n")
    print("$(rank): $(local_to_global(indices))"); print("\n")
end 
Uh=TrialFESpace(Vh)

# Define integration mesh and quadrature
order=1
# Define manufactured functions
u(x) = x[1]+x[2]^order
f(x) = -Δ(u)(x)
degree = 2*order+1
Ω = Triangulation(rdmodel)
dΩ = Measure(Ω,degree)

a(u,v) = ∫( ∇(v)⊙∇(u) )*dΩ
b(v) = ∫(v*f)*dΩ

op = AffineFEOperator(a,b,Uh,Vh)
uh = solve(op)
e = u - uh

# # Compute errors
el2 = sqrt(sum( ∫( e*e )*dΩ ))
eh1 = sqrt(sum( ∫( e*e + ∇(e)⋅∇(e) )*dΩ ))

tol=1e-8
@assert el2 < tol
@assert eh1 < tol

# Uh=FESpace(rdmodel,reffe,conformity=:H1;dirichlet_tags=Int[])
# vector_partition=map(partition(UH.gids)) do indices
#     ones(local_length(indices))
# end
# fv_UH=PVector(vector_partition,partition(UH.gids))
# uH=FEFunction(UH,fv_UH)
# uh=interpolate(uH,Uh)