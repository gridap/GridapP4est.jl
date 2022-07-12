function Base.map(::typeof(Gridap.Arrays.testitem),
  a::Tuple{<:AbstractVector{<:AbstractVector{<:VectorValue}},<:AbstractVector{<:Gridap.Fields.LinearCombinationFieldVector}})
a2=Gridap.Arrays.testitem(a[2])
a1=Vector{eltype(eltype(a[1]))}(undef,size(a2,1))
a1.=zero(Gridap.Arrays.testitem(a1))
(a1,a2)
end

# Fixes Err3 (see below)
function Gridap.Geometry.is_change_possible(
strian::Gridap.Geometry.Triangulation,
ttrian::Gridap.Geometry.Triangulation)
if strian === ttrian || num_cells(strian)==num_cells(ttrian)==0
return true
end
Gridap.Helpers.@check get_background_model(strian) === get_background_model(ttrian) "Triangulations do not point to the same background discrete model!"
D = num_cell_dims(strian)
sglue = get_glue(strian,Val(D))
tglue = get_glue(ttrian,Val(D))
Gridap.Geometry.is_change_possible(sglue,tglue) # Fails here
end

# Fixes Err3 (see below)
function Gridap.CellData.change_domain(a::CellField,
                        ::ReferenceDomain,
                        ttrian::Triangulation,
                        ::ReferenceDomain)
msg = """\n
We cannot move the given CellField to the reference domain of the requested triangulation.
Make sure that the given triangulation is either the same as the triangulation on which the
CellField is defined, or that the latter triangulation is the background of the former.
"""
strian = get_triangulation(a)
if strian === ttrian || num_cells(strian)==num_cells(ttrian)==0
return a
end
@assert Gridap.Geometry.is_change_possible(strian,ttrian) msg
D = num_cell_dims(strian)
sglue = get_glue(strian,Val(D))
tglue = get_glue(ttrian,Val(D))
Gridap.CellData.change_domain_ref_ref(a,ttrian,sglue,tglue)
end

function Gridap.FESpaces.get_cell_fe_data(fun,f,ttrian)
  sface_to_data = fun(f)
  strian = get_triangulation(f)
  if strian === ttrian || num_cells(strian)==num_cells(ttrian)==0
    return sface_to_data
  end
  @assert Gridap.Geometry.is_change_possible(strian,ttrian)
  D = num_cell_dims(strian)
  sglue = get_glue(strian,Val(D))
  tglue = get_glue(ttrian,Val(D))
  Gridap.FESpaces.get_cell_fe_data(fun,sface_to_data,sglue,tglue)
end

function Gridap.Geometry.best_target(trian1::Triangulation,trian2::Triangulation)
  if (num_cells(trian1)==num_cells(trian2)==0)
    return trian1
  end
  @check Gridap.Geometry.is_change_possible(trian1,trian2)
  @check Gridap.Geometry.is_change_possible(trian2,trian1)
  D1 = num_cell_dims(trian1)
  D2 = num_cell_dims(trian2)
  glue1 = get_glue(trian1,Val(D2))
  glue2 = get_glue(trian2,Val(D1))
  Gridap.Geometry.best_target(trian1,trian2,glue1,glue2)
end
