module GridapP4est
  using FillArrays
  using Gridap
  using PartitionedArrays
  using GridapDistributed
  using P4est_wrapper
  const PArrays = PartitionedArrays
  include("UniformlyRefinedForestOfOctreesDiscreteModels.jl")
  include("CubedSphereDistributedDiscreteModels.jl")
  export UniformlyRefinedForestOfOctreesDiscreteModel
end
