
using Gridap, PartitionedArrays, GridapDistributed
using GridapP4est

np = 1
parts = with_mpi() do distribute
  distribute(LinearIndices((prod(np),)))
end

cmodel = CartesianDiscreteModel((0,1,0,1),(1,1))
model  = OctreeDistributedDiscreteModel(parts,cmodel,2) 

reffe = ReferenceFE(raviart_thomas,Float64,1)
Vh = FESpace(model,reffe;dirichlet_tags="boundary") # Fails

cell_dof_ids = map(get_cell_dof_ids,local_views(Vh))

map(parts,local_views(Vh)) do p,Vh
  dof_ids = get_cell_dof_ids(Vh)
  sleep(p*3.0)
  println(">>> Part ", p , ", ndofs=",num_free_dofs(Vh), ", type = ", typeof(Vh))
  display(dof_ids)
end
