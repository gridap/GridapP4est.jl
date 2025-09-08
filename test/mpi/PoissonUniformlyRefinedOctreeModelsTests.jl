
using MPI
using PartitionedArrays

include("../PoissonUniformlyRefinedOctreeModelsTests.jl")
import .PoissonUniformlyRefinedOctreeModelsTests as TestModule

if !MPI.Initialized()
  MPI.Init()
end

parsed_args = TestModule.parse_commandline()
subdomains = Tuple(parsed_args["coarse-cells"])
num_uniform_refinements = parsed_args["num-uniform-refinements"]
num_ghost_layers = parsed_args["num-ghost-layers"]

with_mpi() do distribute
  TestModule.run(distribute, subdomains, num_uniform_refinements, num_ghost_layers)
end

MPI.Finalize()

