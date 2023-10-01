function Base.map(::typeof(Gridap.Arrays.testitem),
  a::Tuple{<:AbstractVector{<:AbstractVector{<:VectorValue}},<:AbstractVector{<:Gridap.Fields.LinearCombinationFieldVector}})
  a2=Gridap.Arrays.testitem(a[2])
  a1=Vector{eltype(eltype(a[1]))}(undef,size(a2,1))
  a1.=zero(Gridap.Arrays.testitem(a1))
  (a1,a2)
end

# Required to transfer fine-grid VECTOR-VALUED fields into coarse-grid 
function Gridap.Adaptivity.FineToCoarseField(fine_fields::AbstractArray{<:Gridap.Fields.Field},
                                              rrule::Gridap.Adaptivity.RefinementRule,
                                              child_ids::AbstractArray{<:Integer})
  
  grid=Gridap.Adaptivity.get_ref_grid(rrule)
  D=num_cell_dims(grid)
  x=zero(Point{D,Float64})
  ffx=lazy_map(evaluate,fine_fields,Fill([x],length(fine_fields)))
  ffx=ffx[1]
  fields = Vector{Gridap.Fields.Field}(undef,Gridap.Adaptivity.num_subcells(rrule))
  fields = fill!(fields,Gridap.Fields.ConstantField(zero(eltype(ffx))))
  for (k,id) in enumerate(child_ids)
    fields[id] = fine_fields[k]
  end
  return Gridap.Adaptivity.FineToCoarseField(fields,rrule)
end