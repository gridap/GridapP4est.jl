using MPI
using PartitionedArrays

include("../issue_108.jl")
import .Issue108Tests as TestModule

if !MPI.Initialized()
  MPI.Init()
end

with_mpi() do distribute
  TestModule.run(distribute)
end

MPI.Finalize()
