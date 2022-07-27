
struct DistributedPatchDecomposition{Dc,Dp,A,B} <: GridapType
  patch_decompositions::A
  model::B
end

function PatchDecomposition(model::DistributedDiscreteModel{Dc,Dp}; Dr=0) where {Dc,Dp}
  patch_decompositions=map_parts(model.models) do model
    PatchDecomposition(model;Dr=Dr)
  end
  A=eltype(patch_decompositions)
  B=eltype(model)
  DistributedPatchDecomposition{Dc,Dp,A,B}(patch_decompositions,model)
end

function Gridap.Geometry.Triangulation(a::DistributedPatchDecomposition)
  trians=map_parts(a.patch_decompositions) do a
    Triangulation(a)
  end
  GridapDistributed.DistributedTriangulation(trians,a.model)
end
