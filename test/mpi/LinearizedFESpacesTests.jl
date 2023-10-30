
using MPI
using PartitionedArrays

include("../LinearizedFESpacesTests.jl")
import .LinearizedFESpacesTests as TestModule

if !MPI.Initialized()
  MPI.Init()
end


with_mpi() do distribute 
  TestModule.run(distribute)
end 

MPI.Finalize()

