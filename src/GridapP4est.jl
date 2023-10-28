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

  include("Environment.jl")
  include("UniformlyRefinedForestOfOctreesDiscreteModels.jl")
  include("OctreeDistributedDiscreteModels.jl")
  include("GridapFixes.jl")
  include("FESpaces.jl")
  include("LinearizedFESpaces.jl")
  include("AdaptivityFlagsMarkingStrategies.jl")
  
  export UniformlyRefinedForestOfOctreesDiscreteModel
  export OctreeDistributedDiscreteModel
  export adapt
  export refine
  export coarsen
  export redistribute
  export nothing_flag, refine_flag, coarsen_flag
  export FixedFractionAdaptiveFlagsMarkingStrategy
  export update_adaptivity_flags!
  
end
