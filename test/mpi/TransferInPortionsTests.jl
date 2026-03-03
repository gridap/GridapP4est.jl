
using MPI
using PartitionedArrays

include("../TransferInPortionsTests.jl")
import .TransferInPortionsTests as TestModule

if !MPI.Initialized()
  MPI.Init()
end

with_mpi() do distribute 
  TestModule.run(distribute)
end 

MPI.Finalize()
