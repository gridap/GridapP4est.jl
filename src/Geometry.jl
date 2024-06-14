struct FaceToCellGlueNonConforming{A,B,C,D} <: GridapType
  face_is_owner    :: Vector{Bool}
  face_to_subface  :: Vector{Int8}
  face_to_cell_glue:: Gridap.Geometry.FaceToCellGlue{A,B,C,D}
end

function _adjust_cell_to_lface_to_pindex!(f2cg_nc,Df,ncglue,cell_faces)
  f2cg=f2cg_nc.face_to_cell_glue 
  face_to_cell=f2cg.face_to_cell
  face_to_lface=f2cg.face_to_lface
  cell_to_lface_to_pindex=f2cg.cell_to_lface_to_pindex
  for i=1:length(f2cg_nc.face_is_owner)
    if (f2cg_nc.face_is_owner[i])
      cell=face_to_cell[i]
      lface=face_to_lface[i]
      gface=cell_faces[cell][lface]
      (oface,_,_)=ncglue.owner_faces_lids[Df][gface]
      pindex=ncglue.owner_faces_pindex[Df][oface]
      s=f2cg.cell_to_lface_to_pindex.ptrs[cell]
      f2cg.cell_to_lface_to_pindex.data[s+lface-1]=pindex
    end 
  end
end 

function FaceToCellGlueNonConforming(topo,
                                     cell_grid,
                                     face_grid,
                                     face_to_bgface,
                                     bgface_to_lcell,
                                     face_is_owner,
                                     face_to_subface,
                                     ncglue)

    face_to_cell_glue=Gridap.Geometry.FaceToCellGlue(topo,
                                     cell_grid,
                                     face_grid,
                                     face_to_bgface,
                                     bgface_to_lcell)

    f2cg_nc=FaceToCellGlueNonConforming(face_is_owner,
                                        face_to_subface,
                                        face_to_cell_glue)

    Dc=num_cell_dims(cell_grid)
    Df=num_cell_dims(face_grid)
    cell_faces=Gridap.Geometry.get_faces(topo,Dc,Df)
    _adjust_cell_to_lface_to_pindex!(f2cg_nc,Df,ncglue,cell_faces)
    f2cg_nc     
end 


function _generate_owner_face_to_subfaces(model,ncglue::GridapP4est.NonConformingGlue{D}) where D
  num_children=GridapP4est.get_num_children(Val{D-1})
  
  # Generate owner_face_to_subfaces_ptrs
  num_owner_faces=length(keys(ncglue.owner_faces_lids[D-1]))
  owner_face_to_subfaces_ptrs  = Vector{Int}(undef, num_owner_faces+1)
  owner_face_to_subfaces_ptrs .= num_children
  length_to_ptrs!(owner_face_to_subfaces_ptrs)
  
  topology=Gridap.Geometry.get_grid_topology(model)
  cell_faces=Gridap.Geometry.get_faces(topology,D,D-1)

  num_hanging_faces = ncglue.num_hanging_faces[D]
  num_regular_faces = ncglue.num_regular_faces[D]

  # Generate owner_face_to_subfaces_data
  owner_face_to_subfaces_data=Vector{Int}(undef,owner_face_to_subfaces_ptrs[end]-1)
  owner_face_to_subfaces_data.=-1
  for i=1:num_hanging_faces
    (ocell,ocell_lface,subface)=ncglue.hanging_faces_glue[D][i]
    if (ocell!=-1)
      ocell_lface_within_dim =GridapP4est.face_lid_within_dim(Val{D}, ocell_lface)
      owner_gface=cell_faces[ocell][ocell_lface_within_dim]
      (lowner,_,_)=ncglue.owner_faces_lids[D-1][owner_gface]
      spos=owner_face_to_subfaces_ptrs[lowner]
      owner_face_to_subfaces_data[spos+subface-1]=num_regular_faces+i
    end 
  end 
  Gridap.Arrays.Table(owner_face_to_subfaces_data,owner_face_to_subfaces_ptrs)
end 

function Gridap.Geometry.SkeletonTriangulation(model::DiscreteModel{D},
                                               ncglue::GridapP4est.NonConformingGlue{D}) where {D}

      num_regular_faces=ncglue.num_regular_faces[D]
      owner_faces_lids=ncglue.owner_faces_lids

      topo=Gridap.Geometry.get_grid_topology(model)
      is_boundary_face=Gridap.Geometry.get_isboundary_face(topo,D-1)

      # Include all regular faces wich are 
      # either owner or not in the boundary
      face_to_regular_bgface=Int[]
      for (gface,ibf) in enumerate(is_boundary_face)
         if gface <= num_regular_faces
            if (!ibf)
              push!(face_to_regular_bgface,gface)
            else 
              if (haskey(owner_faces_lids[D-1],gface))
                push!(face_to_regular_bgface,gface)
              end 
            end 
         end 
      end

      # IMPORTANT NOTE: Plus side contains the hanging side of all cell interfaces, while 
      #                 Minus side the "owner" side of all cell interfaces. This is a MUST.
      #                 For facet integration purposes, the Jacobian of the geometrical mapping
      #                 is extracted from the plus side. 

      # Plus side  
      regular_bgface_to_lcell_minus=Vector{Int8}(undef,num_regular_faces)   
      regular_bgface_to_lcell_minus.=2
      for i=1:length(is_boundary_face)
        if ( i<=num_regular_faces && 
              is_boundary_face[i]  && 
              !(haskey(owner_faces_lids[D-1],i)))
          regular_bgface_to_lcell_minus[i]=1
        end
      end
      plus=BoundaryTriangulation(
        model,
        face_to_regular_bgface,
        regular_bgface_to_lcell_minus,
        ncglue)
      
      
      # Minus side
      regular_bgface_to_lcell_plus=Fill(Int8(1),num_regular_faces)
      minus=BoundaryTriangulation(
        model,
        face_to_regular_bgface,
        regular_bgface_to_lcell_plus,
        ncglue)
  

      SkeletonTriangulation(plus,minus)
end 

function Gridap.Geometry.BoundaryTriangulation(
  model::DiscreteModel{D},
  face_to_regular_bgface::AbstractVector{<:Integer},
  regular_bgface_to_lcell::AbstractVector{<:Integer},
  ncglue::GridapP4est.NonConformingGlue{D}) where D

  @assert length(regular_bgface_to_lcell)==ncglue.num_regular_faces[D]

  T=eltype(face_to_regular_bgface)
  face_to_bgface=T[]
  face_is_owner=Bool[]
  face_to_subface=Int8[]
  bgface_to_lcell=Vector{T}(undef,
                            ncglue.num_regular_faces[D]+ncglue.num_hanging_faces[D])

  bgface_to_lcell.=1                    

  num_children=GridapP4est.get_num_children(Val{D-1})
  owner_face_to_subfaces=_generate_owner_face_to_subfaces(model,ncglue)

  @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))]: owner_face_to_subfaces=$(owner_face_to_subfaces)"

  # Generate face_to_bgface AND bgface_to_lcell 
  for (i,gface) in enumerate(face_to_regular_bgface)
    if (haskey(ncglue.owner_faces_lids[D-1],gface))
      # Find all subfaces of current owner face!
      (lowner,_,_)=ncglue.owner_faces_lids[D-1][gface]
      subfaces=owner_face_to_subfaces[lowner]

      # Regular face is owner of several other hanging faces
      if (regular_bgface_to_lcell[gface]==1)
        bgface_to_lcell[gface]=1
        for i=1:num_children
          if (subfaces[i]!=-1)
            push!(face_to_bgface, gface)    
            push!(face_is_owner,true)
            push!(face_to_subface,i)
          end 
        end
      else
        for i=1:num_children
          if (subfaces[i]!=-1)
            push!(face_to_bgface, subfaces[i])
            bgface_to_lcell[subfaces[i]]=1
            push!(face_is_owner,false)
            push!(face_to_subface,-1)
          end
        end
      end 
    else
      # Regular face which is not owner of any other hanging face
      push!(face_to_bgface, gface)
      bgface_to_lcell[gface]=regular_bgface_to_lcell[gface]
      push!(face_is_owner,false)
      push!(face_to_subface,-1)
    end 
  end 

  topo = Gridap.Geometry.get_grid_topology(model)
  cell_grid = get_grid(model)
  bgface_grid = Gridap.Geometry.Grid(ReferenceFE{D-1},model)
  bgface_grid = view(bgface_grid,face_to_bgface)

  glue=FaceToCellGlueNonConforming(topo,
                              cell_grid,
                              bgface_grid,
                              face_to_bgface,
                              bgface_to_lcell,
                              face_is_owner,
                              face_to_subface,
                              ncglue)

  trian = Gridap.Geometry.BodyFittedTriangulation(model,bgface_grid,face_to_bgface)
  BoundaryTriangulation(trian,glue)
end

function Gridap.Geometry.SkeletonTriangulation(model::Gridap.Adaptivity.AdaptedDiscreteModel{D},
                                               ncglue::GridapP4est.NonConformingGlue{D}) where {D}
  trian = SkeletonTriangulation(Gridap.Adaptivity.get_model(model),ncglue)
  return Gridap.Adaptivity.AdaptedTriangulation(trian,model)
end

function Gridap.Geometry.SkeletonTriangulation(
  portion,model::OctreeDistributedDiscreteModel{Dc};kwargs...) where Dc
  gids   = get_face_gids(model.dmodel,Dc)
  trians = map(local_views(model.dmodel),
               partition(gids),
               model.non_conforming_glue) do model, gids, ncglue
    trian=SkeletonTriangulation(model,ncglue)
    GridapDistributed.filter_cells_when_needed(portion,gids,trian)
  end
  GridapDistributed.DistributedTriangulation(trians,model)
end

struct SubfaceRefCoordsMap{A} <: Gridap.Arrays.Map
  face_is_owner    :: Vector{Bool}
  face_to_subface  :: Vector{Int8}
  ref_rule         :: A
end 

function Gridap.Arrays.return_cache(f::SubfaceRefCoordsMap,gface::Integer)
  ref_grid=Gridap.Adaptivity.get_ref_grid(f.ref_rule)
  poly=Gridap.Adaptivity.get_polytope(f.ref_rule)
  poly_vertices=Gridap.Geometry.get_vertex_coordinates(poly)
  cell_coordinates=get_cell_coordinates(ref_grid)
  cache_cell_coordinates=array_cache(cell_coordinates)
  poly_vertices,cell_coordinates,cache_cell_coordinates
end

function Gridap.Arrays.evaluate!(cache,f::SubfaceRefCoordsMap,gface::Integer)
  poly_vertices,cell_coordinates,cache_cell_coordinates=cache
  if (f.face_is_owner[gface])
    subface=f.face_to_subface[gface]
    getindex!(cache_cell_coordinates,cell_coordinates,subface)
  else 
    poly_vertices
  end   
end 

function Gridap.Geometry.get_glue(trian::BoundaryTriangulation{Dc,Dp,A,<:FaceToCellGlueNonConforming},
                                  ::Val{D},::Val{D}) where {D,Dc,Dp,A}
  tface_to_mface = trian.glue.face_to_cell_glue.face_to_cell
  
  face_to_q_vertex_coords = 
      Gridap.Geometry._compute_face_to_q_vertex_coords(trian,trian.glue.face_to_cell_glue)
  f(p) = 
  Gridap.ReferenceFEs.get_shapefuns(Gridap.ReferenceFEs.LagrangianRefFE(Float64,Gridap.ReferenceFEs.get_polytope(p),1))
  ftype_to_shapefuns = map( f, Gridap.Geometry.get_reffes(trian) )
  face_to_shapefuns = Gridap.ReferenceFEs.expand_cell_data(ftype_to_shapefuns,trian.glue.face_to_cell_glue.face_to_ftype)
  face_s_q = lazy_map(Gridap.Fields.linear_combination,face_to_q_vertex_coords,face_to_shapefuns)

  # Map subface to cell ref space
  ref_rule=nothing
  if Dc==1
    ref_rule=Gridap.Adaptivity.RefinementRule(SEGMENT,2)
  else 
    @assert Dc==2
    ref_rule=Gridap.Adaptivity.RedRefinementRule(QUAD)
  end 

  m=SubfaceRefCoordsMap(trian.glue.face_is_owner,trian.glue.face_to_subface,ref_rule)
  subface_ref_coords = lazy_map(m, Gridap.Arrays.IdentityVector(length(trian.glue.face_to_subface)) )
  face_to_q_vertex_coords = lazy_map(evaluate,face_s_q,subface_ref_coords) 
  face_s_q = lazy_map(Gridap.Fields.linear_combination,face_to_q_vertex_coords,face_to_shapefuns)
    
  tface_to_mface_map = face_s_q
  mface_to_tface = nothing
  Gridap.Geometry.FaceToFaceGlue(tface_to_mface,tface_to_mface_map,mface_to_tface)
end

function Gridap.Geometry.get_facet_normal(trian::BoundaryTriangulation{Dc,Dp,A,<:FaceToCellGlueNonConforming}) where {Dc,Dp,A}
  Gridap.Geometry.get_facet_normal(trian, trian.glue.face_to_cell_glue)
end 
