
using Test
using Gridap, Gridap.Geometry, Gridap.Adaptivity, Gridap.ReferenceFEs
using PartitionedArrays, GridapDistributed
using GridapP4est

function same_model(m1::GridapDistributed.DistributedDiscreteModel,m2)
  res = map(local_views(m1),local_views(m2)) do m1,m2
    same_model(m1,m2)
  end
  return reduce(&,res)
end

function same_model(m1::DiscreteModel{Dc},m2) where Dc
  t1 = get_grid_topology(m1)
  t2 = get_grid_topology(m2)

  mcoords1 = get_node_coordinates(m1)
  mcoords2 = get_node_coordinates(m2)
  mids1 = get_cell_node_ids(m1)
  mids2 = get_cell_node_ids(m2)
  mcellcoords1 = lazy_map(Broadcasting(Reindex(mcoords1)),mids1)
  mcellcoords2 = lazy_map(Broadcasting(Reindex(mcoords2)),mids2)
  A = all(map(x -> x in mcoords2, mcoords1))
  
  tcoords1 = Geometry.get_vertex_coordinates(t1)
  tcoords2 = Geometry.get_vertex_coordinates(t2)
  tids1 = Geometry.get_faces(t1,Dc,0)
  tids2 = Geometry.get_faces(t2,Dc,0)
  tcellcoords1 = lazy_map(Reindex(tcoords1),tids1.data)
  tcellcoords2 = lazy_map(Reindex(tcoords2),tids2.data)
  B = all(map(x -> x in tcoords2, tcoords1))

  return A && B
end

function run(ranks,Dc,isperiodic)
  nc = (Dc==2) ? (3,3) : (3,3,3)
  domain = (Dc==2) ? (0,3,0,3) : (0,3,0,3,0,3)
  good_cmodel = CartesianDiscreteModel(domain,nc;isperiodic=isperiodic)
  d_cmodel = OctreeDistributedDiscreteModel(ranks,good_cmodel,0)
  
  cmodel = PartitionedArrays.getany(local_views(d_cmodel)) 
  @test same_model(good_cmodel,cmodel)
  
  d_fmodel, d_glue = refine(d_cmodel);
  fmodel = PartitionedArrays.getany(local_views(d_fmodel))
  
  good_fmodel = CartesianDiscreteModel(domain,nc.*2;isperiodic=isperiodic)
  @test same_model(good_fmodel,fmodel)
end

ranks = with_mpi() do distribute
  distribute(LinearIndices((1,)))
end

for isperiodic in [(true,true),(true,false),(false,true)]
  run(ranks,2,isperiodic)
end

for isperiodic in [
  (true,true,true),
  (false,true,true),(true,false,true),(true,true,false),
  (false,false,true),(true,false,false),(false,true,false)
  ]
  run(ranks,3,isperiodic)
end
