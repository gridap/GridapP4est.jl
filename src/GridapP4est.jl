module GridapP4est
  using FillArrays
  using Gridap
  using Gridap.FESpaces
  using Gridap.Geometry
  using Gridap.Adaptivity
  using Gridap.ReferenceFEs

  using MPI
  using PartitionedArrays
  using GridapDistributed
  using P4est_wrapper
  const PArrays = PartitionedArrays

  include("PartitionedArraysExtensions.jl")
  include("Environment.jl")
  include("UniformlyRefinedForestOfOctreesDiscreteModels.jl")
  include("OctreeDistributedDiscreteModels.jl")
  
  export UniformlyRefinedForestOfOctreesDiscreteModel
  export OctreeDistributedDiscreteModel
  export refine
  export coarsen
  export redistribute
  export octree_distributed_discrete_model_free!

  export generate_node_coordinates
  export setup_cell_prange
  export generate_grid_and_topology
end
