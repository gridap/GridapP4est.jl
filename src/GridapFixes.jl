function Base.map(::typeof(Gridap.Arrays.testitem),
  a::Tuple{<:AbstractVector{<:AbstractVector{<:VectorValue}},<:AbstractVector{<:Gridap.Fields.LinearCombinationFieldVector}})
  a2=Gridap.Arrays.testitem(a[2])
  a1=Vector{eltype(eltype(a[1]))}(undef,size(a2,1))
  a1.=zero(Gridap.Arrays.testitem(a1))
  (a1,a2)
end

"""
  Given a AdaptivityGlue and a CellField defined on the parent(old) mesh, 
  returns an equivalent CellField on the child(new) mesh.
"""
function Gridap.Adaptivity.change_domain_o2n(
     f_old,
     old_trian::Triangulation{Dc},
     new_trian::AdaptedTriangulation,
     glue::AdaptivityGlue{<:Gridap.Adaptivity.MixedGlue}) where {Dc}

  oglue = get_glue(old_trian,Val(Dc))
  nglue = get_glue(new_trian,Val(Dc))
    
  Gridap.Helpers.@notimplementedif num_point_dims(old_trian) != Dc
  Gridap.Helpers.@notimplementedif isa(nglue,Nothing)

  if (num_cells(old_trian) != 0)
    n2o_faces_map=glue.n2o_faces_map[end]
    new_coarsened_cells_ids=findall(i->n2o_faces_map.ptrs[i+1]-n2o_faces_map.ptrs[i]>1,
                                    1:length(n2o_faces_map.ptrs)-1)
    # If mixed refinement/coarsening, then f_c2f is a Table
    f_old_data=Gridap.CellData.get_data(f_old)
    f_c2f=Gridap.Adaptivity.c2f_reindex(f_old_data,glue)
    new_rrules = Gridap.Adaptivity.get_new_cell_refinement_rules(glue)
    field_array=lazy_map(OldToNewField, f_c2f, new_rrules, glue.n2o_cell_to_child_id)
    return Gridap.CellData.similar_cell_field(f_old,field_array,new_trian,ReferenceDomain())
  else
    f_new = Fill(Gridap.Fields.ConstantField(0.0),num_cells(new_trian))
    return CellData.similar_cell_field(f_old,f_new,new_trian,ReferenceDomain())
  end 
end 

function Gridap.Adaptivity.get_new_cell_refinement_rules(g::AdaptivityGlue{<:Gridap.Adaptivity.MixedGlue})
  old_rrules = g.refinement_rules
  n2o_faces_map = g.n2o_faces_map[end]
  @assert isa(n2o_faces_map,Gridap.Arrays.Table)
  lazy_map(Reindex(old_rrules), [n2o_faces_map.data[n2o_faces_map.ptrs[i]] for i=1:length(n2o_faces_map.ptrs)-1])
end


"""
Given a domain and a non-overlapping refined cover, a `FineToCoarseField`
is a `Field` defined in the domain and constructed by a set of fields defined on 
the subparts of the covering partition.
The refined cover is represented by a `RefinementRule`. 
"""

abstract type NewFieldType end;
struct CoarsenedNewFieldType <: NewFieldType end;
struct RefinedOrUntouchedNewFieldType <: NewFieldType end;   

# Unfortunately, I cannot use traits in OldToNewField as its 
# usage results in conversion errors when leveraging it with 
# the lazy_map infrastructure
struct OldToNewField <: Gridap.Fields.Field
  new_field_type::NewFieldType
  fine_to_coarse_field
  refined_or_untouched_field
  function OldToNewField(a::NewFieldType,b::Gridap.Fields.Field,c::Gridap.Fields.Field)
    new(a,b,c)
  end 
end

function OldToNewField(old_fields::AbstractArray{<:Gridap.Fields.Field},
                       rrule::RefinementRule,
                       child_ids::AbstractVector{<:Integer})
  println("child_ids: $(child_ids)")
  @assert length(old_fields)==length(child_ids)                   
  if length(old_fields)==1
    cell_map = get_cell_map(rrule)[child_ids[1]]
    old_field=old_fields[1]
    fine_to_coarse_field=Gridap.Adaptivity.FineToCoarseField(
          [old_field for i=1:Gridap.Adaptivity.num_subcells(rrule)],
          rrule,
          [i for i=1:Gridap.Adaptivity.num_subcells(rrule)])
    refined_or_untouched_field=old_field∘cell_map
    OldToNewField(RefinedOrUntouchedNewFieldType(),fine_to_coarse_field,refined_or_untouched_field)
  else 
    @assert length(old_fields) <= Gridap.Adaptivity.num_subcells(rrule)
    fine_to_coarse_field=Gridap.Adaptivity.FineToCoarseField(old_fields,rrule,child_ids)
    cell_map = get_cell_map(rrule)[1]
    refined_or_untouched_field=old_fields[1]∘cell_map 
    OldToNewField(CoarsenedNewFieldType(),fine_to_coarse_field,refined_or_untouched_field)
  end
end

function Gridap.Fields.return_cache(a::OldToNewField,x::AbstractArray{<:Point})
  f2c_cache=Gridap.Fields.return_cache(a.fine_to_coarse_field,x)
  rou_cache=Gridap.Fields.return_cache(a.refined_or_untouched_field,x)
  return (f2c_cache,rou_cache)
end

function Gridap.Fields.evaluate!(cache,a::OldToNewField,x::AbstractArray{<:Point})
  if isa(a.new_field_type,CoarsenedNewFieldType)
    f2c_cache,rou_cache = cache
    Gridap.Fields.evaluate!(f2c_cache,a.fine_to_coarse_field,x)
  else
    @assert isa(a.new_field_type,RefinedOrUntouchedNewFieldType)
    f2c_cache,rou_cache = cache
    Gridap.Fields.evaluate!(rou_cache,a.refined_or_untouched_field,x)
  end  
end