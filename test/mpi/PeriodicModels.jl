using MPI
using PartitionedArrays

include("../PeriodicModels.jl")
import .PeriodicModelsTests as TestModule

if !MPI.Initialized()
  MPI.Init()
end

with_mpi() do distribute 
  TestModule.run(distribute)
end 

MPI.Finalize()
