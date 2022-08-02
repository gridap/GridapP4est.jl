module GridapP4est
  using FillArrays
  using Gridap
  using Gridap.FESpaces
  using Gridap.Geometry

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

  include("UniformlyRefinedForestOfOctreesDiscreteModels.jl")
  export UniformlyRefinedForestOfOctreesDiscreteModel

  include("GMG/PartitionedArraysExtensions.jl")
  include("GMG/OctreeDistributedDiscreteModels.jl")
  include("GMG/RedistributeTools.jl")
  include("GMG/ModelHierarchies.jl")
  include("GMG/FESpaceHierarchies.jl")
  include("GMG/InterGridTransferOperators.jl")
  include("GMG/GridapFixes.jl")
  include("GMG/RichardsonSmoothers.jl")
  include("GMG/JacobiLinearSolvers.jl")
  include("GMG/GMGLinearSolvers.jl")

  include("GMG/PatchBasedSmoothers/seq/PatchDecompositions.jl")
  include("GMG/PatchBasedSmoothers/seq/PatchFESpaces.jl")
  include("GMG/PatchBasedSmoothers/seq/PatchBasedLinearSolvers.jl")
  include("GMG/PatchBasedSmoothers/mpi/PatchDecompositions.jl")
  include("GMG/PatchBasedSmoothers/mpi/PatchFESpaces.jl")

  export change_domain_fine_to_coarse
  export change_domain_coarse_to_fine

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
  export GMGLinearSolver
  export RichardsonSmoother
  export JacobiLinearSolver
  export PatchDecomposition
  export PatchFESpace
  export PatchBasedLinearSolver

end
