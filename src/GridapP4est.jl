module GridapP4est
  using FillArrays
  using Gridap

  # BEG TO-REMOVE
  using Gridap.Helpers: @check
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
  include("RedistributeTools.jl")
  include("ModelHierarchies.jl")
  include("FESpaceHierarchies.jl")
  include("InterGridTransferOperators.jl")
  include("GridapFixes.jl")
  include("JacobiSmoothers.jl")
  include("GMGLinearSolvers.jl")

  export change_domain_fine_to_coarse
  export change_domain_coarse_to_fine

  export UniformlyRefinedForestOfOctreesDiscreteModel
  export OctreeDistributedDiscreteModel
  export refine
  export redistribute
  export FineToCoarseModelGlue
  export octree_distributed_discrete_model_free!

  # InterGridTransferOperators
  export InterpolationMat
  export RestrictionMat
  export InterpolationMat
  export RestrictionMat
  export setup_interpolations_and_restrictions

  # ModelHierarchy
  export ModelHierarchy
  export ModelHierarchyLevel
  export model_hierarchy_free!
  export num_levels
  export get_level
  export get_level_model
  export get_level_model_before_redist

  # FESpaceHierarchy
  export FESpaceHierarchyLevel
  export redistribute_fe_function
  export get_level_fe_space
  export get_level_fe_space_before_redist

  # Solvers
  export GMG!
  export JacobiSmoother

end
