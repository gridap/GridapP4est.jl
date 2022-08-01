module MeshHierarchiesTests
  using MPI
  using PartitionedArrays
  using Gridap
  using GridapP4est
  using P4est_wrapper

  function run(parts,num_parts_x_level)
    domain=(0,1,0,1)
    partition=(1,1)
    cmodel=CartesianDiscreteModel(domain,partition)
    mh=ModelHierarchy(parts,cmodel,num_parts_x_level)
    model_hierarchy_free!(mh)
  end

  # Give me how many processors you want per level
  # in an array with as many entries as levels
  num_parts_x_level = [4,2,2,2,2,2]
  if !MPI.Initialized()
    MPI.Init()
  end
  parts = get_part_ids(mpi,6)
  run(parts,num_parts_x_level)
  MPI.Finalize()
end
