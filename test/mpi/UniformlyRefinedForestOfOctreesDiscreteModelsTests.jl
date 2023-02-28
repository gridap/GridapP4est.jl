
using MPI
using PartitionedArrays

include("../UniformlyRefinedForestOfOctreesDiscreteModelsTests.jl")
import .UniformlyRefinedForestOfOctreesDiscreteModelsTests as TestModule

if !MPI.Initialized()
  MPI.Init()
end

parsed_args = TestModule.parse_commandline()
subdomains = Tuple(parsed_args["subdomains"])
num_uniform_refinements = parsed_args["num-uniform-refinements"]

with_backend(TestModule.run,MPIBackend(),(prod(subdomains)),subdomains,num_uniform_refinements)

MPI.Finalize()

