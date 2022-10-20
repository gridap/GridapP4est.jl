module GridapP4est
  using FillArrays
  using Gridap
  using Gridap.FESpaces
  using Gridap.Geometry
  using Gridap.Refinement

  # BEG TO-REMOVE
  using MPI
  using LinearAlgebra
  using IterativeSolvers
  using Printf
  # END TO-REMOVE

  using PartitionedArrays
  using GridapDistributed
  using P4est_wrapper
  const PArrays = PartitionedArrays

  include("PartitionedArraysExtensions.jl")
  include("UniformlyRefinedForestOfOctreesDiscreteModels.jl")
  include("OctreeDistributedDiscreteModels.jl")
  
  export UniformlyRefinedForestOfOctreesDiscreteModel
  export OctreeDistributedDiscreteModel
  export refine
  export redistribute
  export octree_distributed_discrete_model_free!
end
