function Gridap.CellData.change_domain(a::CellField,target_trian::Triangulation,target_domain::DomainStyle)
  strian=get_triangulation(a) 
  if (strian===target_trian)
    change_domain(a,DomainStyle(a),target_domain)
  else
    change_domain(a,get_triangulation(a),DomainStyle(a),target_trian,target_domain)
  end
end

# I copied this fix as-is from GridapSolvers.jl !!! We should decide where to put it.
function Gridap.Adaptivity.change_domain_n2o(f_fine,ftrian::Gridap.Adaptivity.AdaptedTriangulation{Dc},ctrian::Gridap.Geometry.Triangulation,glue::Gridap.Adaptivity.AdaptivityGlue{<:Gridap.Adaptivity.RefinementGlue}) where Dc
  fglue = Gridap.Geometry.get_glue(ftrian,Val(Dc))
  cglue = Gridap.Geometry.get_glue(ctrian,Val(Dc))

  Gridap.Helpers.@notimplementedif Gridap.Geometry.num_point_dims(ftrian) != Dc
  Gridap.Helpers.@notimplementedif isa(cglue,Nothing)

  if (num_cells(ctrian) != 0)
    ### New Triangulation -> New Model
    fine_tface_to_field = Gridap.CellData.get_data(f_fine)
    fine_mface_to_field = Gridap.Geometry.extend(fine_tface_to_field,fglue.mface_to_tface)

    ### New Model -> Old Model
    # f_c2f[i_coarse] = [f_fine[i_fine_1], ..., f_fine[i_fine_nChildren]]
    f_c2f = Gridap.Adaptivity.f2c_reindex(fine_mface_to_field,glue)

    child_ids = Gridap.Adaptivity.f2c_reindex(glue.n2o_cell_to_child_id,glue)
    rrules    = Gridap.Adaptivity.get_old_cell_refinement_rules(glue)
    coarse_mface_to_field = lazy_map(Gridap.Adaptivity.FineToCoarseField,f_c2f,rrules,child_ids)

    ### Old Model -> Old Triangulation
    coarse_tface_to_field = lazy_map(Reindex(coarse_mface_to_field),cglue.tface_to_mface)
    f_coarse = lazy_map(Broadcasting(∘),coarse_tface_to_field,cglue.tface_to_mface_map)

    return Gridap.CellData.similar_cell_field(f_fine,f_coarse,ctrian,ReferenceDomain())
  else
    f_coarse = Fill(Gridap.Fields.ConstantField(0.0),num_cells(fcoarse))
    return Gridap.CellData.similar_cell_field(f_fine,f_coarse,ctrian,ReferenceDomain())
  end
end

function Base.map(::typeof(Gridap.Arrays.testitem),
  a::Tuple{<:AbstractVector{<:AbstractVector{<:VectorValue}},<:AbstractVector{<:Gridap.Fields.LinearCombinationFieldVector}})
  a2=Gridap.Arrays.testitem(a[2])
  a1=Vector{eltype(eltype(a[1]))}(undef,size(a2,1))
  a1.=zero(Gridap.Arrays.testitem(a1))
  (a1,a2)
end


"""
  Given a map from new cells to old cells, computes the inverse map.
  In the general case (refinement+coarsening), the n2o map is a Table, 
  but the algorithm is optimized for Vectors (refinement only).

  IMPORTANT NOTE: This can go to Gridap.Adaptivity as far proper tests are written.

"""
function Gridap.Adaptivity.get_o2n_faces_map(ncell_to_ocell::Gridap.Arrays.Table{T}) where {T<:Integer}
  println(ncell_to_ocell.data)
  nC = maximum(ncell_to_ocell.data)
  ptrs = fill(0,nC+1)
  for ccell in ncell_to_ocell.data
    ptrs[ccell+1] += 1
  end
  Gridap.Arrays.length_to_ptrs!(ptrs)
  data = Vector{Int}(undef,ptrs[end]-1)
  for fcell=1:length(ncell_to_ocell.ptrs)-1
    for j=ncell_to_ocell.ptrs[fcell]:ncell_to_ocell.ptrs[fcell+1]-1
      ccell=ncell_to_ocell.data[j]
      data[ptrs[ccell]]=fcell
      ptrs[ccell]+=1
    end
  end
  Gridap.Arrays.rewind_ptrs!(ptrs)
  ocell_to_ncell=Gridap.Arrays.Table(data,ptrs)
  println(ocell_to_ncell.ptrs)
  println(ocell_to_ncell.data)
  ocell_to_ncell
end

"""
  Given a AdaptivityGlue and a CellField defined on the parent(old) mesh, 
  returns an equivalent CellField on the child(new) mesh.
"""
function Gridap.Adaptivity.change_domain_o2n(
     f_old,
     old_trian::Triangulation{Dc},
     new_trian::AdaptedTriangulation,
     glue::AdaptivityGlue) where {Dc}

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

function OldToNewField(old_fields::AbstractArray{<:Gridap.Fields.Field},rrule::RefinementRule,child_id::Integer)
  if length(old_fields)==1
    cell_map = get_cell_map(rrule)[child_id]
    old_field=old_fields[1]
    fine_to_coarse_field=Gridap.Adaptivity.FineToCoarseField(
          [old_field for i=1:Gridap.Adaptivity.num_subcells(rrule)],
          rrule)
    refined_or_untouched_field=old_field∘cell_map
    OldToNewField(RefinedOrUntouchedNewFieldType(),fine_to_coarse_field,refined_or_untouched_field)
  else 
    Gridap.Helpers.@check length(old_fields) == Gridap.Adaptivity.num_subcells(rrule)
    Gridap.Helpers.@check child_id==-1
    cell_map = get_cell_map(rrule)[1]
    fine_to_coarse_field=Gridap.Adaptivity.FineToCoarseField(old_fields,rrule)
    refined_or_untouched_field=old_fields[1]∘cell_map 
    OldToNewField(CoarsenedNewFieldType(),fine_to_coarse_field,refined_or_untouched_field)
  end
end

function OldToNewField(old_field::Gridap.Fields.Field,rrule::RefinementRule,child_id::Integer)
  Gridap.Helpers.@notimplemented
end

# # Necessary for distributed meshes, where not all children of a coarse cell may belong to the processor. 
# function FineToCoarseField(fine_fields::AbstractArray{<:Field},rrule::RefinementRule,child_ids::AbstractArray{<:Integer})
#   fields = Vector{Field}(undef,num_subcells(rrule))
#   fields = fill!(fields,ConstantField(0.0))
#   for (k,id) in enumerate(child_ids)
#     fields[id] = fine_fields[k]
#   end
#   return FineToCoarseField(fields,rrule)
# end

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

# Fast evaluation of FineToCoarseFields: 
# Points are pre-classified into the children cells, which allows for the search to be 
# skipped entirely. 
function Gridap.Fields.return_cache(a::OldToNewField,x::AbstractArray{<:Point},child_ids::AbstractArray{<:Integer})
  Gridap.Helpers.@notimplemented
end

function Gridap.Fields.evaluate!(cache,a::OldToNewField,x::AbstractArray{<:Point},child_ids::AbstractArray{<:Integer})
  Gridap.Helpers.@notimplemented
end