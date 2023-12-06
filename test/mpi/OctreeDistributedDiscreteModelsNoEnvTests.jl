
using MPI
using PartitionedArrays

include("../OctreeDistributedDiscreteModelsTests.jl")
import .OctreeDistributedDiscreteModelsTests as TestModule

MPI.Init()

# 2D 
with_mpi() do distribute
  TestModule.run(distribute,4,(2,2),[4,2])
end

# 3D
with_mpi() do distribute
  TestModule.run(distribute,4,(2,2,2),[4,2])
end

MPI.Finalize()
