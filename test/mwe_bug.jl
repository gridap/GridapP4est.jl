using P4est_wrapper
using GridapP4est
using Gridap
using PartitionedArrays

np = 1
ranks = with_mpi() do distribute
  distribute(LinearIndices((np,)))
end

Dc = 3
domain = (Dc==2) ? (0,1,0,1) : (0,1,0,1,0,1)
nc     = (Dc==2) ? (2,2) : (2,2,2)

nrefs  = 1
cmodel = CartesianDiscreteModel(domain,nc)
model  = UniformlyRefinedForestOfOctreesDiscreteModel(ranks,cmodel,nrefs) # Breaks

############################################################################################
# Going low-level on UniformlyRefinedForestOfOctreesDiscreteModel() constructor

import GridapP4est: setup_ptr_pXest_objects

comm = ranks.comm
ptr_pXest_connectivity,ptr_pXest,ptr_pXest_ghost,ptr_pXest_lnodes = setup_ptr_pXest_objects(Val{Dc},
                                                                                            comm,
                                                                                            cmodel,
                                                                                            nrefs)

pXest = ptr_pXest[]
p8est_tree_array_index(pXest.trees, 0)[] # Fails

Int(pXest.trees[].elem_size) # 216
sizeof(p8est_tree_t) # 168
sizeof(p4est_tree_t) # 216
# So it seems we have a p4est_tree_t instead? 
