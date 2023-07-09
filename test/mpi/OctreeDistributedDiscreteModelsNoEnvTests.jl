
using MPI
using PartitionedArrays

include("../OctreeDistributedDiscreteModelsTests.jl")
import .OctreeDistributedDiscreteModelsTests as TestModule

MPI.Init()

with_mpi() do distribute
  TestModule.run(distribute,(2,2),[4,2])
end 

MPI.Finalize()
