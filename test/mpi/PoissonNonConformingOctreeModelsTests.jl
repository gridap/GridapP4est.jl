
using MPI
using PartitionedArrays

include("../PoissonNonConformingOctreeModelsTests.jl")
import .PoissonNonConformingOctreeModelsTests as TestModule

if !MPI.Initialized()
  MPI.Init()
end

with_mpi() do distribute 
  TestModule.run(distribute)
end 

MPI.Finalize()

