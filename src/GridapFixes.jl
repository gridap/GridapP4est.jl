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