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
Uh=FESpace(rdmodel,reffe,conformity=:H1;dirichlet_tags=Int[])
# Uh=FESpace(rdmodel,reffe,conformity=:H1;dirichlet_tags=Int[])
# vector_partition=map(partition(UH.gids)) do indices
#     ones(local_length(indices))
# end
# fv_UH=PVector(vector_partition,partition(UH.gids))
# uH=FEFunction(UH,fv_UH)
# uh=interpolate(uH,Uh)