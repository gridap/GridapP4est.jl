module CubedSphereDistributedDiscreteModelsTests
  using GridapP4est
  using PartitionedArrays
  using Test
  prun(mpi,4) do parts
    num_uniform_refinements=4
    model=CubedSphereDiscreteModel(parts,num_uniform_refinements)
  end
end
