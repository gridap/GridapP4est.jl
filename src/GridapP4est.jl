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
  include("PXestTypeMethods.jl")
  include("UniformlyRefinedForestOfOctreesDiscreteModels.jl")
  include("OctreeDistributedDiscreteModels.jl")
  include("Geometry.jl")
  include("AnisotropicallyAdapted3DDistributedDiscreteModels.jl")
  include("GridapFixes.jl")
  include("FESpaces.jl")
  include("AdaptivityFlagsMarkingStrategies.jl")

  
  export UniformlyRefinedForestOfOctreesDiscreteModel
  export OctreeDistributedDiscreteModel
  export AnisotropicallyAdapted3DDistributedDiscreteModel
  export vertically_adapt 
  export horizontally_adapt
  export vertically_uniformly_refine
  export nothing_flag, refine_flag, coarsen_flag
  export FixedFractionAdaptiveFlagsMarkingStrategy
  export update_adaptivity_flags!
  
end
