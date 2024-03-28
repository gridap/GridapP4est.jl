
using MPI
using PartitionedArrays

include("../PoissonAnisotropicOctreeModelsTests.jl")
import .PoissonAnisotropicOctreeModelsTests as TestModule

if !MPI.Initialized()
  MPI.Init()
end

with_mpi() do distribute 
  TestModule.run(distribute)
end 

MPI.Finalize()

