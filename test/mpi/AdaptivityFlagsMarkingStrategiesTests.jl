
using MPI
using PartitionedArrays

include("../AdaptivityFlagsMarkingStrategiesTests.jl")
import .AdaptivityFlagsMarkingStrategiesTests as TestModule

if !MPI.Initialized()
  MPI.Init()
end

with_mpi() do distribute 
  TestModule.run(distribute)
end 

MPI.Finalize()

