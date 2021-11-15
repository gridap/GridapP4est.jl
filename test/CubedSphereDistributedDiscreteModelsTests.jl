module CubedSphereDistributedDiscreteModelsTests
  using GridapP4est
  using Gridap
  using GridapDistributed
  using PartitionedArrays
  using Test
  prun(mpi,4) do parts
    num_uniform_refinements=4
    model=CubedSphereDiscreteModel(parts,num_uniform_refinements)
    Ω=Triangulation(model)
    writevtk(Ω,"Ω")
  end
end
