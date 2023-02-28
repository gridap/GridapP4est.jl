
using MPI
using PartitionedArrays

include("../OctreeDistributedDiscreteModelsTests.jl")
import .OctreeDistributedDiscreteModelsTests as TestModule

with_backend(TestModule.run_with_env,MPIBackend(),4,(2,2),[4,2])
MPI.Finalize()
