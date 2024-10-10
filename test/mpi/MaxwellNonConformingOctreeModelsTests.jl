
using MPI
using PartitionedArrays

include("../MaxwellNonConformingOctreeModelsTests.jl")
import .MaxwellNonConformingOctreeModelsTests as TestModule

if !MPI.Initialized()
  MPI.Init()
end

with_mpi() do distribute 
  TestModule.run(distribute)
end 

MPI.Finalize()

