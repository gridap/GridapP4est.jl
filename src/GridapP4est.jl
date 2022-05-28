module GridapP4est
  using FillArrays
  using Gridap
  using PartitionedArrays
  using GridapDistributed
  using P4est_wrapper
  const PArrays = PartitionedArrays
  include("UniformlyRefinedForestOfOctreesDiscreteModels.jl")
  include("OctreeDistributedDiscreteModels.jl")
  export UniformlyRefinedForestOfOctreesDiscreteModel
  export OctreeDistributedDiscreteModel
  export refine
  export octree_distributed_discrete_model_free
end
