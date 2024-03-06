function Base.map(::typeof(Gridap.Arrays.testitem),
  a::Tuple{<:AbstractVector{<:AbstractVector{<:VectorValue}},<:AbstractVector{<:Gridap.Fields.LinearCombinationFieldVector}})
  a2=Gridap.Arrays.testitem(a[2])
  a1=Vector{eltype(eltype(a[1]))}(undef,size(a2,1))
  a1.=zero(Gridap.Arrays.testitem(a1))
  (a1,a2)
end
