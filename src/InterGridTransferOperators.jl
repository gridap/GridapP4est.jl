function change_domain_coarse_to_fine(c_cell_field,
                                      ftrian::Triangulation{Dc,Dp},
                                      glue::Union{Nothing,FineToCoarseModelGlue}) where {Dc,Dp}

  if (num_cells(ftrian) != 0)
    fcell_to_child_id=glue.fcell_to_child_id
    fcell_to_ccell=glue.fine_to_coarse_faces_map[Dc+1]

    Gridap.Helpers.@check length(fcell_to_child_id)==num_cells(ftrian)
    Gridap.Helpers.@check DomainStyle(c_cell_field)==ReferenceDomain()

    # TO-DO: can we pre-compute this?
    rrule_reffe=Gridap.ReferenceFEs.LagrangianRefFE(Float64,QUAD,1)
    rrule_grid=Gridap.Geometry.UnstructuredGrid(
                  Gridap.Visualization.compute_reference_grid(rrule_reffe,2))

    rrule_f_to_c_ref_map=get_cell_map(rrule_grid)
    cfield=Gridap.CellData.get_data(c_cell_field)
    m1=Reindex(cfield)
    f_to_cfield=lazy_map(m1,fcell_to_ccell)
    m2=Reindex(rrule_f_to_c_ref_map)
    cfield_to_ffield=lazy_map(∘,f_to_cfield,lazy_map(m2,fcell_to_child_id))
    Gridap.CellData.GenericCellField(cfield_to_ffield,ftrian,ReferenceDomain())
  else
    cfield_to_ffield=Fill(Gridap.Fields.ConstantField(0.0),num_cells(ftrian))
    Gridap.CellData.GenericCellField(cfield_to_ffield,ftrian,ReferenceDomain())
  end
end

function change_domain_fine_to_coarse(f_cell_field,
                                      ctrian::Triangulation{Dc,Dp},
                                      glue::FineToCoarseModelGlue) where {Dc,Dp}
  Gridap.Helpers.@check DomainStyle(f_cell_field)==ReferenceDomain()
  # BEG TO-DO: can we pre-compute something of the following?
  rrule_reffe=Gridap.ReferenceFEs.LagrangianRefFE(Float64,QUAD,1)
  rrule_grid=Gridap.Geometry.UnstructuredGrid(
                Gridap.Visualization.compute_reference_grid(rrule_reffe,2))
  rrule_f_to_c_ref_map=get_cell_map(rrule_grid)
  rrule_c_to_f_ref_map=lazy_map(Gridap.Fields.inverse_map,rrule_f_to_c_ref_map)
  model=Gridap.Geometry.UnstructuredDiscreteModel(rrule_grid)
  Ω=Triangulation(model)
  cache1 = Gridap.CellData._point_to_cell_cache(Gridap.CellData.KDTreeSearch(),Ω)
  x_to_cell(x) = Gridap.CellData._point_to_cell!(cache1, x)
  ffield=Gridap.CellData.get_data(f_cell_field)
  A=typeof(ffield)
  B=typeof(x_to_cell)
  C=typeof(rrule_c_to_f_ref_map)
  m=TransformFineFieldsToCoarseFieldMap{Dc,A,B,C}(ffield,x_to_cell,rrule_c_to_f_ref_map)
  c_cell_field=lazy_map(m,Gridap.Arrays.IdentityVector(num_cells(ctrian)))
  Gridap.CellData.GenericCellField(c_cell_field,ctrian,ReferenceDomain())
  # END TO-DO: can we pre-compute something of the following?
end

struct TransformFineFieldsToCoarseFieldMap{Dc,A,B,C} <: Gridap.Fields.Map
  ffield               :: A
  x_to_cell            :: B
  rrule_c_to_f_ref_map :: C
end

struct FineFieldsToCoarseField{Dc,A,B,C,D} <: Gridap.Fields.Field
  ccell                :: A
  ffield               :: B
  x_to_cell            :: C
  rrule_c_to_f_ref_map :: D
end


function Gridap.Arrays.return_cache(m::TransformFineFieldsToCoarseFieldMap,ccell::Integer)
  nothing
end

function Gridap.Arrays.evaluate!(cache,
                    m::TransformFineFieldsToCoarseFieldMap{Dc,A,B,C},
                    ccell::Integer) where {Dc,A,B,C}
  FineFieldsToCoarseField{Dc,typeof(ccell),A,B,C}(ccell,
                                                m.ffield,
                                                m.x_to_cell,
                                                m.rrule_c_to_f_ref_map)
end

function Gridap.Arrays.return_cache(m::FineFieldsToCoarseField{Dc},
                      x::AbstractVector{<:Point}) where Dc

  cffield=array_cache(m.ffield)
  cx=array_cache(x)
  xi=getindex!(cx,x,1)
  fi=getindex!(cffield,m.ffield,1)
  ti=Gridap.Arrays.return_type(fi,xi)
  cfx=Gridap.Arrays.CachedArray(ti,1)
  cf=Gridap.Arrays.return_cache(fi,xi)
  crrule_map=array_cache(m.rrule_c_to_f_ref_map)
  crrule_map_i=getindex!(crrule_map,m.rrule_c_to_f_ref_map,1)
  crrule_map_i_cache=Gridap.Arrays.return_cache(crrule_map_i,xi)
  x_to_f=collect(x)
  x_to_child=Vector{UInt8}(undef,length(x))
  for i in eachindex(x)
    xi=getindex!(cx,x,i)
    child_id=m.x_to_cell(xi)
    x_to_child[i]=child_id
    crrule_map_i=getindex!(crrule_map,m.rrule_c_to_f_ref_map,child_id)
    x_to_f[i]=evaluate!(crrule_map_i_cache,crrule_map_i,xi)
  end
  cffield, cx, cf, cfx, x, x_to_f, x_to_child
end

function Gridap.Arrays.evaluate!(cache,
                    m::FineFieldsToCoarseField{Dc},
                    x::AbstractVector{<:Point}) where Dc
  cffield, cx, cf, cfx, xc, x_to_f, x_to_child=cache
  Gridap.Helpers.@notimplementedif !(xc === x)
  Gridap.Arrays.setsize!(cfx, size(x))
  num_children=2^Dc
  base = (m.ccell-1)*num_children
  for i in eachindex(x)
    xi=x_to_f[i]
    f=getindex!(cffield,m.ffield,base+x_to_child[i])
    cfx.array[i]=evaluate!(cf,f,xi)
  end
  cfx.array
end

function change_domain_coarse_to_fine(Uh::GridapDistributed.DistributedFESpace,
                                      c_cell_field,
                                      ftrian::GridapDistributed.DistributedTriangulation{Dc,Dp},
                                      glue::MPIData{<:Union{Nothing,FineToCoarseModelGlue}}) where {Dc,Dp}

  i_am_in_coarse=(c_cell_field != nothing)

  fields=map_parts(GridapDistributed.local_views(ftrian)) do Ω
    if (i_am_in_coarse)
      c_cell_field.fields.part
    else
      Gridap.Helpers.@check num_cells(Ω) == 0
      Gridap.CellData.GenericCellField(Fill(Gridap.Fields.ConstantField(0.0),num_cells(Ω)),Ω,ReferenceDomain())
    end
  end
  c_cell_field_fine=GridapDistributed.DistributedCellField(fields)

  dfield=map_parts(change_domain_coarse_to_fine,
                    GridapDistributed.local_views(c_cell_field_fine),
                    GridapDistributed.local_views(ftrian),
                    glue)
  GridapDistributed.DistributedCellField(dfield)
end

function ProlongationOperator(uH,Uh,Vh,ftrian,glue,dΩ)
  uH_h = change_domain_coarse_to_fine(uH,ftrian,glue)
  a(u,v) = ∫(v*u)dΩ
  l(v)   = ∫(v*uH_h)dΩ
  op=AffineFEOperator(a,l,Uh,Vh)
  solve(op)
end

function change_domain_fine_to_coarse(f_cell_field::GridapDistributed.DistributedCellField,
                                      ctrian::Union{GridapDistributed.DistributedTriangulation,Nothing},
                                      glue)
  i_am_in_coarse=(ctrian != nothing)
  println("??? $(i_am_in_coarse) $(glue.part)" )
  if i_am_in_coarse
    c_f_cell_field,cglue=map_parts(GridapDistributed.local_views(ctrian)) do _
      f_cell_field.fields.part,glue.part
    end

    dfield=map_parts(change_domain_fine_to_coarse,
                     c_f_cell_field,
                     GridapDistributed.local_views(ctrian),
                     cglue)
    GridapDistributed.DistributedCellField(dfield)
  else
    return nothing
  end
end

function RestrictionOperator(uh,UH,VH,ctrian,glue,dΩ)
  uh_H = change_domain_fine_to_coarse(uh,ctrian,glue)
  a(u,v) = ∫(v*u)dΩ
  l(v)   = ∫(v*uh_H)dΩ
  op=AffineFEOperator(a,l,UH,VH)
  solve(op)
end
