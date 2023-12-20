
using MPI
using PartitionedArrays

include("../OctreeDistributedDiscreteModelsTests.jl")
import .OctreeDistributedDiscreteModelsTests as TestModule

# 2D 
with_mpi() do distribute
  TestModule.run_with_env(distribute,4,(2,2),[4,2])
end

# 3D
with_mpi() do distribute
  TestModule.run_with_env(distribute,4,(2,2,2),[4,2])
end

MPI.Finalize()
