
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

# num_ghost_layers preserved across refine/adapt/coarsen/redistribute
with_mpi() do distribute
  TestModule.run_num_ghost_layers_with_env(distribute)
end

MPI.Finalize()
