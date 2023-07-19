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
  include("RedistributeTools.jl")
  
  export UniformlyRefinedForestOfOctreesDiscreteModel
  export OctreeDistributedDiscreteModel
  export refine
  export coarsen
  export redistribute
  export setup_non_conforming_distributed_discrete_model
  export MPIVoidVector
  export i_am_in
  export nothing_flag, refine_flag, coarsen_flag
  export redistribute_fe_function
  
end
