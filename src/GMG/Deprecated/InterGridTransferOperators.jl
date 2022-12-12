function change_domain_coarse_to_fine(c_cell_field,
                                      ftrian::Triangulation{Dc,Dp},
                                      glue::Union{Nothing,FineToCoarseModelGlue}) where {Dc,Dp}
end

function change_domain_fine_to_coarse(f_cell_field,
                                      ctrian::Triangulation{Dc,Dp},
                                      glue::FineToCoarseModelGlue) where {Dc,Dp}
end

function change_domain_coarse_to_fine(c_cell_field,
                                      ftrian::GridapDistributed.DistributedTriangulation{Dc,Dp},
                                      glue::MPIData{<:Union{Nothing,FineToCoarseModelGlue}}) where {Dc,Dp}
  i_am_in_coarse = (c_cell_field != nothing)

  fields = map_parts(GridapDistributed.local_views(ftrian)) do Ω
    if (i_am_in_coarse)
      c_cell_field.fields.part
    else
      Gridap.Helpers.@check num_cells(Ω) == 0
      Gridap.CellData.GenericCellField(Fill(Gridap.Fields.ConstantField(0.0),num_cells(Ω)),Ω,ReferenceDomain())
    end
  end
  c_cell_field_fine = GridapDistributed.DistributedCellField(fields)

  dfield = map_parts(change_domain_coarse_to_fine,
                    GridapDistributed.local_views(c_cell_field_fine),
                    GridapDistributed.local_views(ftrian), glue)
  return GridapDistributed.DistributedCellField(dfield)
end

function change_domain_fine_to_coarse(f_cell_field::GridapDistributed.DistributedCellField,
                                      ctrian::Union{GridapDistributed.DistributedTriangulation,Nothing},
                                      glue)
  i_am_in_coarse = (ctrian != nothing)
  if i_am_in_coarse
    c_f_cell_field, cglue = map_parts(GridapDistributed.local_views(ctrian)) do _
                                f_cell_field.fields.part, glue.part
                            end

    dfield = map_parts(change_domain_fine_to_coarse,
                     c_f_cell_field,
                     GridapDistributed.local_views(ctrian),
                     cglue)
    return GridapDistributed.DistributedCellField(dfield)
  end
  return nothing
end


mutable struct InterpolationMat
  Ωh
  Ωh_ghost
  dΩh
  UH
  Uh
  Vh
  Uh_red
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

function InterpolationMat(mh,fespaces,level,qdegree)
  Gridap.Helpers.@check 1 <= level <= num_levels(mh)-1

  tests,trials=fespaces

  # Fine level (all processes participate)
  model_h=get_level_model_before_redist(mh,level)
  Ωh=Triangulation(model_h.dmodel)
  Ωh_ghost=Triangulation(with_ghost,model_h.dmodel)
  dΩh=Measure(Ωh,qdegree)


  Vh=get_level_fe_space_before_redist(tests[level])
  Uh=get_level_fe_space_before_redist(trials[level])
  Mh=assemble_mass_matrix(model_h,Vh,Uh,qdegree)


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
    VH=get_level_fe_space(tests[level+1])
    UH=get_level_fe_space(trials[level+1])
    dof_values_H_fe_space_layout = PVector(0.0,VH.gids)
    uH_zero_dirichlet_values=
       map_parts(GridapDistributed.local_views(UH.spaces)) do space
           zeros(num_dirichlet_dofs(space))
       end
  else
    UH = nothing
    dof_values_H_fe_space_layout = nothing
    uH_zero_dirichlet_values = nothing
  end

  cache=InterpolationMat(Ωh, Ωh_ghost, dΩh, UH, Uh, Vh,
                         get_level_fe_space(trials[level]),
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

function LinearAlgebra.mul!(x::PVector,
                            A::InterpolationMat,
                            y::Union{PVector,Nothing};
                            verbose=false,
                            reltol=1.0e-14)
  parts = A.mesh_hierarchy_level.model.parts
  Gridap.Helpers.@check GridapP4est.i_am_in(parts)

  model_next=get_level_model(A.mesh_hierarchy_level_next)
  coarse_parts  = model_next.parts
  i_am_in_coarse= GridapP4est.i_am_in(coarse_parts)

  if (i_am_in_coarse)
    copy!(A.dof_values_H_fe_space_layout,y)
    uH = FEFunction(A.UH,
                    A.dof_values_H_fe_space_layout,
                    A.uH_zero_dirichlet_values)
  else
    uH = nothing
  end

  uH_h = change_domain_coarse_to_fine(uH,A.Ωh_ghost,A.mesh_hierarchy_level.ref_glue)

  l(v) = ∫(v⋅uH_h)A.dΩh
  Gridap.FESpaces.assemble_vector!(l,A.dof_values_h_sys_layout_b,A.Vh)
  parts=A.mesh_hierarchy_level.model.parts

  if (A.mesh_hierarchy_level.model_red != nothing)
    fill!(A.dof_values_h_sys_layout_x,0.0)
    IterativeSolvers.cg!(A.dof_values_h_sys_layout_x,
                         A.mh,
                         A.dof_values_h_sys_layout_b;
                         verbose=(i_am_main(parts) && verbose),
                         reltol=reltol)
    copy!(A.dof_values_h_fe_space_layout,
          A.dof_values_h_sys_layout_x)

    uhold=FEFunction(A.Uh,
                     A.dof_values_h_fe_space_layout,
                     A.uh_zero_dirichlet_values)

    uhnew  = redistribute_fe_function(uhold,
                                      A.Uh,
                                      A.Uh_red,
                                      A.mesh_hierarchy_level.model_red.dmodel,
                                      A.mesh_hierarchy_level.red_glue)

    copy!(x,uhnew.metadata.free_values)
  else
    IterativeSolvers.cg!(x,
                         A.mh,
                         A.dof_values_h_sys_layout_b;
                         verbose=(i_am_main(parts) && verbose),
                         reltol=reltol)
  end
end

struct RestrictionMat
  ΩH
  dΩH
  dΩh
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

function RestrictionMat(mh,fespaces,level,qdegree)
  Gridap.Helpers.@check 1 <= level <= num_levels(mh)-1
  tests,trials=fespaces

  model_h=get_level_model_before_redist(mh,level)
  Ωh  = Triangulation(model_h.dmodel)
  dΩh = Measure(Ωh,qdegree)
  Uh  = get_level_fe_space_before_redist(trials[level])
  Vh  = get_level_fe_space_before_redist(tests[level])
  Vh_red= get_level_fe_space(tests[level])
  Uh_red= get_level_fe_space(trials[level])
  dof_values_h_fe_space_layout_red = PVector(0.0,Vh_red.gids)
  dof_values_h_fe_space_layout     = PVector(0.0,Vh.gids)
  uh_zero_dirichlet_values_red=
    map_parts(GridapDistributed.local_views(Uh_red.spaces)) do space
         zeros(num_dirichlet_dofs(space))
    end

  model_H=get_level_model(mh,level+1)
  if (GridapP4est.i_am_in(model_H.parts))
    ΩH  = Triangulation(model_H.dmodel)
    dΩH = Measure(ΩH,qdegree)
    VH=get_level_fe_space(tests[level+1])
    UH=get_level_fe_space(trials[level+1])
    MH=assemble_mass_matrix(model_H,VH,UH,qdegree)
    dof_values_H_fe_space_layout = PVector(0.0,VH.gids)
    dof_values_H_sys_layout_b = similar(dof_values_H_fe_space_layout,(axes(MH)[1],))
    uH_zero_dirichlet_values=
    map_parts(GridapDistributed.local_views(UH.spaces)) do space
         zeros(num_dirichlet_dofs(space))
    end
  else
    ΩH = nothing
    dΩH = nothing
    VH = nothing
    UH = nothing
    MH = nothing
    dof_values_H_fe_space_layout = nothing
    dof_values_H_sys_layout_b = nothing
    uH_zero_dirichlet_values = nothing
  end

  RestrictionMat(ΩH, dΩH, dΩh, Uh, Uh_red, UH, VH,
                 get_level(mh,level),
                 get_level(mh,level+1),
                 uh_zero_dirichlet_values_red,
                 dof_values_h_fe_space_layout_red,
                 uH_zero_dirichlet_values,
                 dof_values_H_sys_layout_b,
                 MH)
end

function LinearAlgebra.mul!(x::Union{PVector,Nothing},
                            A::RestrictionMat,
                            y::PVector;
                            verbose=false,
                            reltol=1.0e-14)

  parts = A.mesh_hierarchy_level.model.parts
  Gridap.Helpers.@check GridapP4est.i_am_in(parts)

  copy!(A.dof_values_h_fe_space_layout_red,y)
  if (A.mesh_hierarchy_level.model_red != nothing)
    uhold = FEFunction(A.Uh_red,
                       A.dof_values_h_fe_space_layout_red,
                       A.uh_zero_dirichlet_values_red)
                                                                                                                     
    uh=redistribute_fe_function(uhold,
                                A.Uh_red,
                                A.Uh,
                                A.mesh_hierarchy_level.model.dmodel,
                                A.mesh_hierarchy_level.red_glue;
                                reverse=true)
  else
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
    l(v) = ∫(v⋅uh_H)A.dΩH
    Gridap.FESpaces.assemble_vector!(l,A.dof_values_H_sys_layout_b,A.VH)
    fill!(x,0.0)
    IterativeSolvers.cg!(x,
                         A.mH,
                         A.dof_values_H_sys_layout_b;
                         reltol=reltol,
                         verbose=(i_am_main(parts) && verbose))
  end
end

function setup_interpolations_and_restrictions(mh,fespaces,qdegree)
  nlevs          = num_levels(mh)
  interpolations = Vector{InterpolationMat}(undef,nlevs-1)
  restrictions   = Vector{RestrictionMat}(undef,nlevs-1)
  for l = 1:nlevs-1
    model = get_level_model(mh,l)
    if (GridapP4est.i_am_in(model.parts))
      interpolations[l] = InterpolationMat(mh,fespaces,l,qdegree)
      restrictions[l]   = RestrictionMat(mh,fespaces,l,qdegree)
    end
  end
  return interpolations, restrictions
end
