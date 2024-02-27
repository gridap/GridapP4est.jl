
using FillArrays
using Gridap, Gridap.Geometry, Gridap.Adaptivity, Gridap.ReferenceFEs
using PartitionedArrays, GridapDistributed
using GridapP4est

ranks = with_mpi() do distribute
  distribute(LinearIndices((1,)))
end

isperiodic = (true,false)
Dc = 2
good_cmodel = CartesianDiscreteModel((0,3,0,3),(3,3);isperiodic=isperiodic)
d_cmodel = OctreeDistributedDiscreteModel(ranks,good_cmodel,0)

cmodel = PartitionedArrays.getany(local_views(d_cmodel)) 

ctopo = get_grid_topology(cmodel)
good_ctopo = get_grid_topology(good_cmodel)
ctopo_coords = Geometry.get_vertex_coordinates(ctopo)
good_ctopo_coords = Geometry.get_vertex_coordinates(good_ctopo)

d_fmodel, d_glue = refine(d_cmodel);

glue = PartitionedArrays.getany(d_glue)
fmodel = PartitionedArrays.getany(local_views(d_fmodel))

ftopo = get_grid_topology(fmodel)
ftopo_ids = Geometry.get_faces(ftopo,Dc,0)
ftopo_coords = Geometry.get_vertex_coordinates(ftopo)

fmodel_ids = get_cell_node_ids(fmodel)
fmodel_coords = get_node_coordinates(fmodel)

good_fmodel = UnstructuredDiscreteModel(CartesianDiscreteModel((0,3,0,3),(6,6);isperiodic=isperiodic))
good_ftopo = get_grid_topology(good_fmodel)
good_ftopo_ids = Geometry.get_faces(good_ftopo,Dc,0)
good_ftopo_coords = Geometry.get_vertex_coordinates(good_ftopo)
good_fmodel_ids = get_cell_node_ids(good_fmodel)
good_fmodel_coords = get_node_coordinates(good_fmodel)


ftopo_cell_coords = lazy_map(Broadcasting(Reindex(ftopo_coords)),ftopo_ids)
good_ftopo_cell_coords = lazy_map(Broadcasting(Reindex(good_ftopo_coords)),good_ftopo_ids)
