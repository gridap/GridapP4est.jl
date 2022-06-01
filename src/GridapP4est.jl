module GridapP4est
  using FillArrays
  using Gridap

  # BEG TO-REMOVE
  using GridapPETSc
  using Libdl
  using Gridap.Helpers: @check
  using GridapPETSc.PETSC: PetscErrorCode
  using MPI
  # END TO-REMOVE

  using PartitionedArrays
  using GridapDistributed
  using P4est_wrapper
  const PArrays = PartitionedArrays
  include("UniformlyRefinedForestOfOctreesDiscreteModels.jl")
  include("OctreeDistributedDiscreteModels.jl")
  include("InterGridTransferOperators.jl")
  include("PETSCextensions.jl")

  export change_domain_fine_to_coarse
  export change_domain_coarse_to_fine

  export UniformlyRefinedForestOfOctreesDiscreteModel
  export OctreeDistributedDiscreteModel
  export refine
  export FineToCoarseModelGlue
  export octree_distributed_discrete_model_free
  export ProlongationOperator
  export RestrictionOperator

end
