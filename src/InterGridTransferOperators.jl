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

function change_domain_coarse_to_fine(c_cell_field,
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

function change_domain_fine_to_coarse(f_cell_field::GridapDistributed.DistributedCellField,
                                      ctrian::Union{GridapDistributed.DistributedTriangulation,Nothing},
                                      glue)
  i_am_in_coarse=(ctrian != nothing)
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

function assemble_mass_matrix(model,Vh,Uh,degree)
  Ω  = Triangulation(model.dmodel)
  dΩ = Measure(Ω,degree)
  a(u,v)=∫(v⋅u)dΩ
  assemble_matrix(a,Uh,Vh)
end

mutable struct InterpolationMat
  Ωh_red
  dΩh_red
  Ωh
  Ωh_ghost
  ΩH
  dΩh
  UH
  VH
  Uh
  Vh
  Uh_red
  Vh_red
  mesh_hierarchy_level
  mesh_hierarchy_level_next
  uH_zero_dirichlet_values
  dof_values_H_fe_space_layout
  uh_zero_dirichlet_values
  dof_values_h_fe_space_layout
  dof_values_h_sys_layout_b
  dof_values_h_sys_layout_x
  mh
end

# TO-DO: InterpolationMatCache
function LinearAlgebra.mul!(x::PVector,A::InterpolationMat,y::Union{PVector,Nothing})
  parts = A.mesh_hierarchy_level.model.parts
  Gridap.Helpers.@check GridapP4est.i_am_in(parts)

  model_next=get_level_model(A.mesh_hierarchy_level_next)
  coarse_parts  = model_next.parts
  i_am_in_coarse= GridapP4est.i_am_in(coarse_parts)

  if (i_am_in_coarse)
    # To-think: Do we need communication after copy!?
    copy!(A.dof_values_H_fe_space_layout,y)
    exchange!(A.dof_values_H_fe_space_layout)
    uH = FEFunction(A.UH,
                    A.dof_values_H_fe_space_layout,
                    A.uH_zero_dirichlet_values)
    # writevtk(A.ΩH,"uH",nsubcells=2,cellfields=["uH"=>uH])
  else
    uH = nothing
  end

  uH_h = change_domain_coarse_to_fine(uH,A.Ωh_ghost,A.mesh_hierarchy_level.ref_glue)
  #writevtk(A.Ωh,"uH_h",cellfields=["uH_h"=>uH_h])

  l(v) = ∫(v*uH_h)A.dΩh
  Gridap.FESpaces.assemble_vector!(l,A.dof_values_h_sys_layout_b,A.Vh)
  parts=A.mesh_hierarchy_level.model.parts

  if (A.mesh_hierarchy_level.model_red != nothing)
    fill!(A.dof_values_h_sys_layout_x,0.0)
    # println("aaaaaaaaaa aaaaaaaaaaaaaa")
    IterativeSolvers.cg!(A.dof_values_h_sys_layout_x,
                         A.mh,
                         A.dof_values_h_sys_layout_b;
                         verbose=i_am_main(parts),
                         reltol=1.0e-14)
    # println("bbbbbbbbbb bbbbbbbbbbbbbb")
    copy!(A.dof_values_h_fe_space_layout,
          A.dof_values_h_sys_layout_x)
    exchange!(A.dof_values_h_fe_space_layout)

    # map_parts(A.dof_values_h_fe_space_layout.values,
    #           A.dof_values_h_sys_layout_x.values) do a,b
    #   println(a)
    #   println(b)
    # end

    uhold=FEFunction(A.Uh,
                     A.dof_values_h_fe_space_layout,
                     A.uh_zero_dirichlet_values)

    # println("XXX $(sum(∫(uhold)A.dΩh))")

    # map_parts(A.mesh_hierarchy_level.model.parts,uhold.fields) do part,fs
    #     println("AAA1 $(fs.cell_dof_values)")
    #     println("AAA1 $(fs.free_values)")
    #     println("AAA1 $(fs.dirichlet_values)")
    #     println("AAA1 $(fs.fe_space)")
    # end

    # writevtk(A.Ωh,"Correction",cellfields=["correction"=>uhold])

    uhnew  = redistribute_fe_function(uhold,
                                      A.Uh,
                                      A.Uh_red,
                                      A.mesh_hierarchy_level.model_red.dmodel,
                                      A.mesh_hierarchy_level.red_glue)

    #writevtk(A.Ωh_red,"Correction red",cellfields=["correction red"=>uhnew])

    # map_parts(A.mesh_hierarchy_level.model_red.parts,uhnew.fields) do part,fs
    #   if part == 1
    #     println("AAA $(fs.cell_dof_values)")
    #     println("AAA $(fs.free_values)")
    #     println("AAA $(fs.dirichlet_values)")
    #     println("AAA $(fs.fe_space)")
    #   end
    # end

    # println("YYY $(sum(∫(uhnew)A.dΩh_red))")


    # map_parts(uhnew.metadata.free_values.values) do fv
    #   println("RRR $(fv)")
    # end



    # COMM after copy???
    copy!(x,uhnew.metadata.free_values)
    exchange!(x)
    # map_parts(x.values,uhnew.metadata.free_values.values) do x1, x2
    #   println("XXX $(x1)")
    #   println("YYY $(x2)")
    # end
  else
    IterativeSolvers.cg!(x,
                         A.mh,
                         A.dof_values_h_sys_layout_b;
                         verbose=i_am_main(parts),
                         reltol=1.0e-14)
  end
end

function setup_interpolation_mat(mh,fespaces,level)
  Gridap.Helpers.@check 1 <= level <= num_levels(mh)-1

  tests,trials=fespaces
  order=1

  # Fine level (all processes participate)
  model_h=get_level_model_before_redist(mh,level)
  Ωh=Triangulation(model_h.dmodel)
  Ωh_red=Triangulation(get_level_model(mh,level).dmodel)
  Ωh_ghost=Triangulation(with_ghost,model_h.dmodel)
  dΩh=Measure(Ωh,2*(order+1))
  dΩh_red=Measure(Ωh_red,2*(order+1))


  Vh=get_level_fe_space_before_redist(tests[level])
  Uh=get_level_fe_space_before_redist(trials[level])
  Mh=assemble_mass_matrix(model_h,Vh,Uh,2*order+1)


  dof_values_h_fe_space_layout = PVector(0.0,Vh.gids)
  dof_values_h_sys_layout_b = similar(dof_values_h_fe_space_layout,(axes(Mh)[1],))
  dof_values_h_sys_layout_x = similar(dof_values_h_fe_space_layout,(axes(Mh)[2],))
  uh_zero_dirichlet_values=
     map_parts(GridapDistributed.local_views(Uh.spaces)) do space
        zeros(num_dirichlet_dofs(space))
     end

  # Coarse level (only processes of the next level participate)
  model_H=get_level_model(mh,level+1)
  if (GridapP4est.i_am_in(model_H.parts))
    ΩH = Triangulation(model_H.dmodel)
    VH=get_level_fe_space(tests[level+1])
    UH=get_level_fe_space(trials[level+1])
    dof_values_H_fe_space_layout = PVector(0.0,VH.gids)
    uH_zero_dirichlet_values=
       map_parts(GridapDistributed.local_views(UH.spaces)) do space
           zeros(num_dirichlet_dofs(space))
       end
    # map_parts(mh.levels[level].model.parts) do part_id
    #     println("INTERP $(part_id) $(ΩH)")
    # end
  else
    ΩH = nothing
    UH = nothing
    VH = nothing
    dof_values_H_fe_space_layout = nothing
    uH_zero_dirichlet_values = nothing
    # map_parts(mh.levels[level].model.parts) do part_id
    #   println("INTERP $(part_id) $(ΩH)")
    # end
  end

  cache=InterpolationMat(Ωh_red,
                         dΩh_red,
                         Ωh,
                         Ωh_ghost,
                         ΩH,
                         dΩh,
                         UH,
                         VH,
                         Uh,
                         Vh,
                         get_level_fe_space(trials[level]),
                         get_level_fe_space(tests[level]),
                         get_level(mh,level),
                         get_level(mh,level+1),
                         uH_zero_dirichlet_values,
                         dof_values_H_fe_space_layout,
                         uh_zero_dirichlet_values,
                         dof_values_h_fe_space_layout,
                         dof_values_h_sys_layout_b,
                         dof_values_h_sys_layout_x,
                         Mh)
end

mutable struct RestrictionMat
  Ωh
  ΩH
  dΩh
  dΩH
  Uh
  Uh_red
  UH
  VH
  mesh_hierarchy_level
  mesh_hierarchy_level_next
  uh_zero_dirichlet_values_red
  dof_values_h_fe_space_layout_red
  uH_zero_dirichlet_values
  dof_values_H_sys_layout_b
  mH
end

function LinearAlgebra.mul!(x::Union{PVector,Nothing}, A::RestrictionMat, y::PVector)
  parts = A.mesh_hierarchy_level.model.parts
  Gridap.Helpers.@check GridapP4est.i_am_in(parts)

  copy!(A.dof_values_h_fe_space_layout_red,y)
  exchange!(A.dof_values_h_fe_space_layout_red)
  if (A.mesh_hierarchy_level.model_red != nothing)
    uhold = FEFunction(A.Uh_red,
                       A.dof_values_h_fe_space_layout_red,
                       A.uh_zero_dirichlet_values_red)

    # x is in a subcommunicator of y
    uh=redistribute_fe_function(uhold,
                                A.Uh_red,
                                A.Uh,
                                A.mesh_hierarchy_level.model.dmodel,
                                A.mesh_hierarchy_level.red_glue;
                                reverse=true)

    # writevtk(A.Ωh,"Residual",cellfields=["residual"=>uh])
  else
    # x and y are in the same communicator
    uh = FEFunction(A.Uh,
                    A.dof_values_h_fe_space_layout_red,
                    A.uh_zero_dirichlet_values_red)
  end

  # uh is in h communicator, but with void parts for those tasks not in the next level
  # uh_H is void for those tasks not in the next level
  uh_H = change_domain_fine_to_coarse(uh,A.ΩH,A.mesh_hierarchy_level.ref_glue)

  parts = A.mesh_hierarchy_level_next.model.parts
  if (GridapP4est.i_am_in(parts))
    Gridap.Helpers.@check x != nothing
    l(v) = ∫(v*uh_H)A.dΩH
    Gridap.FESpaces.assemble_vector!(l,A.dof_values_H_sys_layout_b,A.VH)
    fill!(x,0.0)
    IterativeSolvers.cg!(x,
                         A.mH,
                         A.dof_values_H_sys_layout_b;
                         verbose=i_am_main(parts),
                         reltol=1.0e-14)
    uH = FEFunction(A.UH,
                    x,
                    A.uH_zero_dirichlet_values)
    # writevtk(A.ΩH,"Residual_projected",cellfields=["residual_projected"=>uH])
  end
end

function setup_restriction_mat(mh,fespaces,level)
  Gridap.Helpers.@check 1 <= level <= num_levels(mh)-1
  tests,trials=fespaces
  order=1

  model_h=get_level_model_before_redist(mh,level)
  Ωh  = Triangulation(model_h.dmodel)
  dΩh = Measure(Ωh,2*order+1)
  Uh  = get_level_fe_space_before_redist(trials[level])
  Vh  = get_level_fe_space_before_redist(tests[level])
  Uh_red= get_level_fe_space(trials[level])
  Vh_red= get_level_fe_space(tests[level])
  dof_values_h_fe_space_layout_red = PVector(0.0,Vh_red.gids)
  dof_values_h_fe_space_layout     = PVector(0.0,Vh.gids)
  uh_zero_dirichlet_values_red=
    map_parts(GridapDistributed.local_views(Uh_red.spaces)) do space
         zeros(num_dirichlet_dofs(space))
    end

  model_H=get_level_model(mh,level+1)
  if (GridapP4est.i_am_in(model_H.parts))
    ΩH  = Triangulation(model_H.dmodel)
    dΩH = Measure(ΩH,2*(order+1))
    VH=get_level_fe_space(tests[level+1])
    UH=get_level_fe_space(trials[level+1])
    MH=assemble_mass_matrix(model_H,VH,UH,2*order+1)
    # Generate arrays of DoF Values for coarse system
    dof_values_H_fe_space_layout = PVector(0.0,VH.gids)
    dof_values_H_sys_layout_b = similar(dof_values_H_fe_space_layout,(axes(MH)[1],))
    uH_zero_dirichlet_values=
    map_parts(GridapDistributed.local_views(UH.spaces)) do space
         zeros(num_dirichlet_dofs(space))
    end
    # map_parts(mh.levels[level].model.parts) do part_id
    #   println("RESTRICT $(part_id) $(ΩH)")
    # end
  else
    ΩH  = nothing
    dΩH = nothing
    VH  = nothing
    UH  = nothing
    MH  = nothing
    dof_values_H_fe_space_layout = nothing
    dof_values_H_sys_layout_b    = nothing
    uH_zero_dirichlet_values = nothing
    # map_parts(mh.levels[level].model.parts) do part_id
    #   println("RESTRICT $(part_id) $(ΩH)")
    # end
  end

  cache=RestrictionMat(Ωh,
                       ΩH,
                       dΩh,
                       dΩH,
                       Uh,
                       Uh_red,
                       UH,
                       VH,
                       get_level(mh,level),
                       get_level(mh,level+1),
                       uh_zero_dirichlet_values_red,
                       dof_values_h_fe_space_layout_red,
                       uH_zero_dirichlet_values,
                       dof_values_H_sys_layout_b,
                       MH)
end
