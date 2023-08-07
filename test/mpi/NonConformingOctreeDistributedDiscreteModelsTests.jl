
using MPI
using PartitionedArrays

include("../NonConformingOctreeDistributedDiscreteModelsTests.jl")
import .NonConformingOctreeDistributedDiscreteModelsTests as TestModule

if !MPI.Initialized()
  MPI.Init()
end

with_mpi() do distribute 
  TestModule.run(distribute)
end 

MPI.Finalize()

