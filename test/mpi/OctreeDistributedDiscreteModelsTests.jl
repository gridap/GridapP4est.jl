
using MPI
using PartitionedArrays

include("../OctreeDistributedDiscreteModelsTests.jl")
import .OctreeDistributedDiscreteModelsTests as TestModule

with_mpi() do distribute
    TestModule.run_with_env(distribute,(2,2),[4,2])
  end 

MPI.Finalize()
