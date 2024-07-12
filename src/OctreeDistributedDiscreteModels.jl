
const nothing_flag  = Cint(0)
const refine_flag   = Cint(1)
const coarsen_flag  = Cint(2)

struct NonConformingGlue{Dc,A,B,C}
  num_regular_faces  :: A # <:AbstractVector{<:AbstractVector{<:Integer}}
  num_hanging_faces  :: B # <:AbstractVector{<:AbstractVector{<:Integer}}
  hanging_faces_glue :: C # <:AbstractVector{<:AbstractVector{<:Tuple{<:Integer,:Integer,Integer}}}}
  function NonConformingGlue(num_regular_faces,num_hanging_faces,hanging_faces_glue)
    Dc = length(num_regular_faces)
    @assert length(num_hanging_faces)==Dc
    @assert length(hanging_faces_glue)==Dc
    A = typeof(num_regular_faces)
    B = typeof(num_hanging_faces)
    C = typeof(hanging_faces_glue)
    new{Dc,A,B,C}(num_regular_faces,num_hanging_faces,hanging_faces_glue)
  end
end


function _create_conforming_model_non_conforming_glue(model::GridapDistributed.DistributedDiscreteModel{Dc}) where Dc
  non_conforming_glue = map(local_views(model)) do model
    num_regular_faces=Vector{Int}(undef,Dc)
    num_hanging_faces=Vector{Int}(undef,Dc)
    hanging_faces_glue=Vector{Tuple{Int,Int,Int}}(undef,Dc)
    for d=1:Dc 
      num_regular_faces[d]=num_faces(model,d-1)
      num_hanging_faces[d]=0
    end
    NonConformingGlue(num_regular_faces,num_hanging_faces,hanging_faces_glue)
  end
  return non_conforming_glue
end

mutable struct OctreeDistributedDiscreteModel{Dc,Dp,A,B,C,D,E,F} <: GridapDistributed.DistributedDiscreteModel{Dc,Dp}
  parts                       :: A
  dmodel                      :: B
  non_conforming_glue         :: C 
  coarse_model                :: D
  ptr_pXest_connectivity      :: E
  ptr_pXest                   :: F
  pXest_type                  :: PXestType
  pXest_refinement_rule_type  :: Union{Nothing,PXestRefinementRuleType}

  # The model for which this variable is true, is the one
  # ultimately responsible for deallocating the pXest_connectivity
  # info
  owns_ptr_pXest_connectivity :: Bool

  # Might be optionally be used, e.g., to enforce that this
  # model is GCed after another existing model
  gc_ref                      :: Any
  function OctreeDistributedDiscreteModel(
    Dc::Int,
    Dp::Int,
    parts,
    dmodel::Union{GridapDistributed.DistributedDiscreteModel,Nothing},
    non_conforming_glue::Union{AbstractVector{<:NonConformingGlue},Nothing},
    coarse_model,
    ptr_pXest_connectivity,
    ptr_pXest,
    pXest_type::PXestType,
    pXest_refinement_rule_type::Union{Nothing,PXestRefinementRuleType},
    owns_ptr_pXest_connectivity::Bool,
    gc_ref)

    if (isa(dmodel,GridapDistributed.DistributedDiscreteModel))
      Gridap.Helpers.@check Dc == Gridap.Geometry.num_cell_dims(dmodel)
      Gridap.Helpers.@check Dp == Gridap.Geometry.num_point_dims(dmodel)
    end

    A = typeof(parts)
    B = typeof(dmodel)
    C = typeof(non_conforming_glue)
    D = typeof(coarse_model)
    E = typeof(ptr_pXest_connectivity)
    F = typeof(ptr_pXest)
    model = new{Dc,Dp,A,B,C,D,E,F}(parts,
                                   dmodel,
                                   non_conforming_glue,
                                   coarse_model,
                                   ptr_pXest_connectivity,
                                   ptr_pXest,
                                   pXest_type,
                                   pXest_refinement_rule_type,
                                   owns_ptr_pXest_connectivity,
                                   gc_ref)
    Init(model)
    return model
  end
end

function OctreeDistributedDiscreteModel(
  parts,
  dmodel::GridapDistributed.DistributedDiscreteModel{Dc,Dp},
  non_conforming_glue::AbstractVector{<:NonConformingGlue{Dc}},
  coarse_model,
  ptr_pXest_connectivity,
  ptr_pXest,
  pXest_type,
  pXest_refinement_rule_type,
  owns_ptr_pXest_connectivity,
  gc_ref) where {Dc,Dp}

  return OctreeDistributedDiscreteModel(Dc,
                                        Dp,
                                        parts,
                                        dmodel,
                                        non_conforming_glue,
                                        coarse_model,
                                        ptr_pXest_connectivity,
                                        ptr_pXest,
                                        pXest_type,
                                        pXest_refinement_rule_type,
                                        owns_ptr_pXest_connectivity,
                                        gc_ref)
end

function _dim_to_pXest_type(Dc)
  Dc==2 ? P4estType() : P8estType() 
end 


function OctreeDistributedDiscreteModel(parts::AbstractVector{<:Integer},
                                        coarse_model::DiscreteModel{Dc,Dp},
                                        num_uniform_refinements) where {Dc,Dp}
  comm = parts.comm
  if GridapDistributed.i_am_in(comm)

    pXest_type = _dim_to_pXest_type(Dc)

    ptr_pXest_connectivity,
      ptr_pXest,
        ptr_pXest_ghost,
          ptr_pXest_lnodes = setup_ptr_pXest_objects(pXest_type,
                                                     comm,
                                                     coarse_model,
                                                     num_uniform_refinements)
    dmodel = setup_distributed_discrete_model(pXest_type,
                                            parts,
                                            coarse_model,
                                            ptr_pXest_connectivity,
                                            ptr_pXest,
                                            ptr_pXest_ghost,
                                            ptr_pXest_lnodes)
    pXest_lnodes_destroy(pXest_type,ptr_pXest_lnodes)
    pXest_ghost_destroy(pXest_type,ptr_pXest_ghost)

    non_conforming_glue = _create_conforming_model_non_conforming_glue(dmodel)

    return OctreeDistributedDiscreteModel(Dc,
                                          Dp,
                                          parts,
                                          dmodel,
                                          non_conforming_glue,
                                          coarse_model,
                                          ptr_pXest_connectivity,
                                          ptr_pXest,
                                          pXest_type,
                                          nothing,
                                          true,
                                          nothing)
  else
    ## HUGE WARNING: Shouldn't we provide here the complementary of parts
    ##               instead of parts? Otherwise, when calling _free!(...)
    ##               we cannot trust on parts.
    return VoidOctreeDistributedDiscreteModel(coarse_model,parts)
  end
end

function OctreeDistributedDiscreteModel(
    parts::AbstractVector{<:Integer},
    coarse_model::DiscreteModel{Dc,Dp}) where {Dc,Dp}
  OctreeDistributedDiscreteModel(parts,coarse_model,0)
end

# Void models

const VoidOctreeDistributedDiscreteModel{Dc,Dp,A,C,D} = OctreeDistributedDiscreteModel{Dc,Dp,A,Nothing,C,D,Nothing}

function VoidOctreeDistributedDiscreteModel(coarse_model::DiscreteModel{Dc,Dp},parts) where {Dc,Dp}
  ptr_pXest_connectivity = setup_pXest_connectivity(coarse_model)
  OctreeDistributedDiscreteModel(Dc,
                                 Dp,
                                 parts,
                                 nothing,
                                 nothing,
                                 coarse_model,
                                 ptr_pXest_connectivity,
                                 nothing,
                                 _dim_to_pXest_type(Dc),
                                 nothing,
                                 true,
                                 nothing)
end

function VoidOctreeDistributedDiscreteModel(model::OctreeDistributedDiscreteModel{Dc,Dp},parts) where {Dc,Dp}
  OctreeDistributedDiscreteModel(Dc,
                                 Dp,
                                 parts,
                                 nothing,
                                 nothing,
                                 model.coarse_model,
                                 model.ptr_pXest_connectivity,
                                 nothing,
                                 _dim_to_pXest_type(Dc),
                                 nothing,
                                 false,
                                 model)
end

# DistributedDiscreteModel API implementation

GridapDistributed.get_parts(model::OctreeDistributedDiscreteModel) = model.parts
GridapDistributed.local_views(model::OctreeDistributedDiscreteModel) = GridapDistributed.local_views(model.dmodel)
GridapDistributed.get_cell_gids(model::OctreeDistributedDiscreteModel) = GridapDistributed.get_cell_gids(model.dmodel)
GridapDistributed.get_face_gids(model::OctreeDistributedDiscreteModel,dim::Integer) = GridapDistributed.get_face_gids(model.dmodel,dim)

# Garbage collection

function octree_distributed_discrete_model_free!(model::VoidOctreeDistributedDiscreteModel{Dc}) where Dc
  if (model.owns_ptr_pXest_connectivity)
    pXest_connectivity_destroy(model.pXest_type,model.ptr_pXest_connectivity)
  end
  return nothing
end

function octree_distributed_discrete_model_free!(model::OctreeDistributedDiscreteModel{Dc}) where Dc
  if !isa(model.ptr_pXest,Nothing)
    pXest_destroy(model.pXest_type,model.ptr_pXest)
  end
  if (model.owns_ptr_pXest_connectivity)
    pXest_connectivity_destroy(model.pXest_type,model.ptr_pXest_connectivity)
  end
  return nothing
end

function Init(a::OctreeDistributedDiscreteModel)
  finalizer(Finalize,a)
end

function Finalize(a::OctreeDistributedDiscreteModel)
  octree_distributed_discrete_model_free!(a)
  return nothing
end

###################################################################
# Private methods

function get_num_children(::P4estType,::PXestUniformRefinementRuleType)
  4
end

function get_num_children(::P8estType,::PXestUniformRefinementRuleType)
  8
end

function get_num_children(::P6estType,::PXestVerticalRefinementRuleType)
  2
end

function get_num_children(::P6estType,::PXestHorizontalRefinementRuleType)
  4
end

function get_num_children(::Type{Val{Dc}}) where Dc
  2^Dc
end

function get_refinement_rule(::P4estType,::PXestUniformRefinementRuleType)
  reffe  = LagrangianRefFE(Float64,QUAD,1)
  Gridap.Adaptivity.RefinementRule(reffe,2)
end 

function get_refinement_rule(::P8estType,::PXestUniformRefinementRuleType)
  reffe  = LagrangianRefFE(Float64,HEX,1)
  Gridap.Adaptivity.RefinementRule(reffe,2)
end 

function get_refinement_rule(::P6estType,::PXestVerticalRefinementRuleType)
  reffe  = LagrangianRefFE(Float64,HEX,1)
  Gridap.Adaptivity.RefinementRule(reffe,(2,1,1))
end 

function get_refinement_rule(::P6estType,::PXestHorizontalRefinementRuleType)
  reffe  = LagrangianRefFE(Float64,HEX,1)
  Gridap.Adaptivity.RefinementRule(reffe,(1,2,2))
end 

function pXest_update_flags!(pXest_type::P4P8estType, ptr_pXest_old, ptr_pXest_new)
  pXest_old = ptr_pXest_old[]
  flags=unsafe_wrap(Array, 
                    Ptr{Cint}(pXest_old.user_pointer), 
                    pXest_old.local_num_quadrants)
  pXest_update_flags!(pXest_type, flags, ptr_pXest_old, ptr_pXest_new)
end

function pXest_update_flags!(pXest_type::P4P8estType, flags, ptr_pXest_old, ptr_pXest_new)
  pXest_old = ptr_pXest_old[]
  pXest_new = ptr_pXest_new[]
  
  num_trees = Cint(pXest_old.connectivity[].num_trees)
  @assert num_trees == Cint(pXest_new.connectivity[].num_trees)

  Dc=num_cell_dims(pXest_type)
  num_children = get_num_children(pXest_type, PXestUniformRefinementRuleType())
  global_iquad_old = 0
  for itree = 0:num_trees-1
   tree_old = pXest_tree_array_index(pXest_type,pXest_old,itree)[]
   tree_new = pXest_tree_array_index(pXest_type,pXest_new,itree)[]
   num_quads_old = Cint(tree_old.quadrants.elem_count)
   local_iquad_old=0
   local_iquad_new=0
   while local_iquad_old < num_quads_old
     q_old = pXest_quadrant_array_index(pXest_type,tree_old,local_iquad_old)
     q_new = pXest_quadrant_array_index(pXest_type,tree_new,local_iquad_new)
     if (pXest_quadrant_is_equal(pXest_type,q_old,q_new))     # q_old was not refined nor coarsened
       flags[global_iquad_old+1] = nothing_flag
       global_iquad_old += 1
       local_iquad_new += 1
       local_iquad_old += 1
     elseif (pXest_quadrant_is_parent(pXest_type,q_old,q_new)) # q_old was refined
       flags[global_iquad_old+1] = refine_flag
       global_iquad_old += 1
       local_iquad_new += num_children
       local_iquad_old += 1
     elseif (pXest_quadrant_is_parent(pXest_type,q_new,q_old)) # q_old and its siblings were coarsened 
       for i=0:num_children-1
         flags[global_iquad_old+i+1] = coarsen_flag
       end
       global_iquad_old += num_children
       local_iquad_old += num_children
       local_iquad_new += 1
     else
       @assert false
     end
   end
  end
end

function p6est_horizontally_adapt_update_flags!(ptr_pXest_old, ptr_pXest_new)
  pXest_old = ptr_pXest_old[]
  pXest_new = ptr_pXest_new[]
  ptr_p4est_old = pXest_old.columns
  ptr_p4est_new = pXest_new.columns
  flags=unsafe_wrap(Array, 
                    Ptr{Cint}(pXest_old.user_pointer), 
                    pXest_old.columns[].local_num_quadrants)
  pXest_update_flags!(P4estType(), flags, ptr_p4est_old, ptr_p4est_new)
end

function p6est_vertically_adapt_update_flags!(ptr_pXest_old, ptr_pXest_new)
  pXest_old = ptr_pXest_old[]
  pXest_new = ptr_pXest_new[]
  flags=unsafe_wrap(Array, 
                    Ptr{Cint}(pXest_old.user_pointer), 
                    pXest_old.layers[].elem_count)
  
  num_trees = Cint(pXest_old.columns[].connectivity[].num_trees)
  @assert num_trees == Cint(pXest_new.columns[].connectivity[].num_trees)

  pXest_type=P6estType()
  num_children = get_num_children(pXest_type, PXestVerticalRefinementRuleType())
  global_icell_old = 0
  for itree = 0:num_trees-1
   tree_old = pXest_tree_array_index(pXest_type,pXest_old,itree)[]
   tree_new = pXest_tree_array_index(pXest_type,pXest_new,itree)[]
   num_quads = Cint(tree_old.quadrants.elem_count)
   num_quads_new = Cint(tree_new.quadrants.elem_count)
   @assert num_quads == num_quads_new

   local_iquad=0
   # Loop over quadrants
   for local_iquad=0:num_quads-1
     q_old = pXest_quadrant_array_index(pXest_type,tree_old,local_iquad)
     q_new = pXest_quadrant_array_index(pXest_type,tree_new,local_iquad)

     f_old,l_old=P6EST_COLUMN_GET_RANGE(q_old[])
     f_new,l_new=P6EST_COLUMN_GET_RANGE(q_new[])
     i_old=f_old
     i_new=f_new
     # Loop over layers within current column
     while i_old<l_old
       q2_old = p2est_quadrant_array_index(pXest_old.layers[], i_old)
       q2_new = p2est_quadrant_array_index(pXest_new.layers[], i_new)
       if (p2est_quadrant_is_equal(q2_old,q2_new))
         flags[global_icell_old+1] = nothing_flag
         i_new += 1
         i_old += 1
         global_icell_old+=1
       elseif (p2est_quadrant_is_ancestor(q2_old,q2_new)) # q2_old was refined
         flags[global_icell_old+1] = refine_flag
         i_new += num_children
         i_old += 1
         global_icell_old+=1
       else # q2_old was coarsened
         @assert p2est_quadrant_is_ancestor(q2_new,q2_old) 
         for i=0:num_children-1
          flags[global_icell_old+i+1] = coarsen_flag
         end
         global_icell_old += num_children
         i_old += num_children
         i_new += 1
       end
     end
    end
  end 
end


function _compute_fine_to_coarse_model_glue(
         cparts,
         cmodel::Union{Nothing,GridapDistributed.DistributedDiscreteModel{Dc}},
         fmodel::GridapDistributed.DistributedDiscreteModel{Dc}) where Dc

  # Fill data for owned (from coarse cells owned by the processor)
  fgids = get_cell_gids(fmodel)
  f1,f2,f3, cgids_snd, cgids_rcv = map(local_views(fmodel),partition(fgids)) do fmodel, fpartition
    if (!(GridapDistributed.i_am_in(cparts)))
      # cmodel might be distributed among less processes than fmodel
      nothing, nothing, nothing, Int[], Int[]
    else
      cgids        = get_cell_gids(cmodel)
      cmodel_local = PArrays.getany(cmodel.models)
      cpartition   = PArrays.getany(partition(cgids))
      fine_to_coarse_faces_map,
        fine_to_coarse_faces_dim,
          fcell_to_child_id =
          _process_owned_cells_fine_to_coarse_model_glue(cmodel_local,fmodel,cpartition,fpartition)
      
      # Note: Reversing snd and rcv    
      lids_rcv,lids_snd = map(PArrays.getany,assembly_local_indices(partition(fgids)))
      @debug "$(MPI.Comm_rank(MPI.COMM_WORLD)): lids_rcv: $(lids_rcv)"
      @debug "$(MPI.Comm_rank(MPI.COMM_WORLD)): lids_snd: $(lids_snd)"
      cgids_data  = local_to_global(cpartition)[fine_to_coarse_faces_map[Dc+1][lids_snd.data]]
      cgids_snd   = PArrays.JaggedArray(cgids_data,lids_snd.ptrs)
      cgids_rcv   = PArrays.JaggedArray(Vector{Int}(undef,length(lids_rcv.data)),lids_rcv.ptrs)

      fine_to_coarse_faces_map, fine_to_coarse_faces_dim, fcell_to_child_id, cgids_snd, cgids_rcv
    end
  end |> tuple_of_arrays 

  # Nearest Neighbors comm: Get data for ghosts (from coarse cells owned by neighboring processors)
  dfcell_to_child_id = map(f3) do fcell_to_child_id
   !isa(fcell_to_child_id,Nothing) ? fcell_to_child_id : Int[]
  end
  cache = fetch_vector_ghost_values_cache(dfcell_to_child_id,partition(fgids))
  fetch_vector_ghost_values!(dfcell_to_child_id,cache) |> wait
  
  # Note: Reversing snd and rcv 
  parts_rcv, parts_snd = assembly_neighbors(partition(fgids))
  PArrays.exchange_fetch!(cgids_rcv,cgids_snd,ExchangeGraph(parts_snd,parts_rcv))

  map(f1,cgids_rcv) do fine_to_coarse_faces_map, cgids_rcv
    if (GridapDistributed.i_am_in(cparts))
      cgids = get_cell_gids(cmodel)
      cpartition = PArrays.getany(partition(cgids))
      # Note: Reversing snd and rcv
      lids_rcv,_ = map(PArrays.getany,assembly_local_indices(partition(fgids)))
      glo_to_loc = global_to_local(cpartition)
      for i in 1:length(lids_rcv.data)
        lid = lids_rcv.data[i]
        gid = cgids_rcv.data[i]
        fine_to_coarse_faces_map[Dc+1][lid] = glo_to_loc[gid]
      end
    end
  end

  # Create distributed glue
  glue = map(f1,f2,f3) do fine_to_coarse_faces_map, fine_to_coarse_faces_dim, fcell_to_child_id
    if (!(GridapDistributed.i_am_in(cparts)))
      nothing
    else
      polytope = (Dc==2 ? QUAD : HEX)
      cmodel_local = PArrays.getany(cmodel.models)
      num_cells_coarse = num_cells(cmodel_local)
      reffe  = LagrangianRefFE(Float64,polytope,1)
      rrules = Fill(Gridap.Adaptivity.RefinementRule(reffe,2),num_cells_coarse)
      AdaptivityGlue(fine_to_coarse_faces_map,fcell_to_child_id,rrules)
    end
  end
  return glue
end

function _process_owned_cells_fine_to_coarse_model_glue(cmodel::DiscreteModel{Dc},
                                                        fmodel::DiscreteModel{Dc},
                                                        cpartition,
                                                        fpartition) where Dc
  fine_to_coarse_faces_map = Vector{Vector{Int}}(undef,Dc+1)
  fine_to_coarse_faces_dim = Vector{Vector{Int}}(undef,Dc)

  num_f_cells   = num_cells(fmodel)                # Number of fine cells (owned+ghost)
  num_o_c_cells = length(own_to_local(cpartition)) # Number of coarse cells (owned)

  # Allocate local vector size # local cells
  fine_to_coarse_faces_map[Dc+1] = Vector{Int}(undef,num_f_cells)
  fcell_to_child_id = Vector{Int}(undef,num_f_cells)

  # Go over all cells of coarse grid portion
  num_children = get_num_children(Val{Dc})
  c = 1
  for cell = 1:num_o_c_cells
    for child = 1:num_children
      fine_to_coarse_faces_map[Dc+1][c+child-1] = cell
      fcell_to_child_id[c+child-1] = child
    end
    c = c + num_children
  end
  ftopology = Gridap.Geometry.get_grid_topology(fmodel)
  for d=1:Dc
    fine_to_coarse_faces_map[d] = Vector{Int}(undef,Gridap.Geometry.num_faces(ftopology,d-1))
    fine_to_coarse_faces_dim[d] = Vector{Int}(undef,Gridap.Geometry.num_faces(ftopology,d-1))
  end
  # c_cell_faces=[]
  # cache_c_cell_faces=[]
  # f_cell_faces=[]
  # cache_f_cell_faces=[]
  # for d=1:Dc
  #   push!(c_cell_faces, Gridap.Geometry.get_faces(ctopology, Dc, d-1))
  #   push!(cache_c_cell_faces,array_cache(last(c_cell_faces)))
  #   push!(f_cell_faces, Gridap.Geometry.get_faces(ftopology, Dc, d-1))
  #   push!(cache_f_cell_faces,array_cache(last(f_cell_faces)))
  # end
  # parent_cell_faces=Vector{Vector{Int}}(undef,Dc)
  # for cell=1:num_f_cells
  #   parent_cell=fine_to_coarse_faces_map[Dc+1][cell]
  #   child=fcell_to_child_id[cell]
  #   for d=1:Dc
  #     parent_cell_faces[d]=getindex!(cache_c_cell_faces[d],
  #                                   c_cell_faces[d],
  #                                   parent_cell)
  #   end
  #   for d=1:Dc
  #     cell_f_faces=getindex!(cache_f_cell_faces[d],
  #                             f_cell_faces[d],
  #                             cell)
  #     for (lf,f) in enumerate(cell_f_faces)
  #       c     = rrule_f_to_c_lid_2D[d][child][lf]
  #       dim_c = rrule_f_to_c_dim_2D[d][child][lf]
  #       if (dim_c == Dc)
  #         fine_to_coarse_faces_map[d][f]=parent_cell
  #       else
  #         fine_to_coarse_faces_map[d][f]=parent_cell_faces[dim_c+1][c]
  #       end
  #       fine_to_coarse_faces_dim[d][f]=dim_c
  #     end
  #   end
  # end

  return fine_to_coarse_faces_map, fine_to_coarse_faces_dim, fcell_to_child_id
end

function _move_fwd_and_check_if_all_children_coarsened(flags,num_o_c_cells,cell,num_children,stride)
  e=cell+num_children*stride-1
  while (cell <= num_o_c_cells) && (cell <= e)
    if (flags[cell]!=coarsen_flag)
      break
    end
    cell = cell+stride
  end
  return cell,cell==e+1
end

function _compute_fine_to_coarse_model_glue(
         pXest_type,
         pXest_refinement_rule_type,
         cparts,
         cmodel::Union{Nothing,GridapDistributed.DistributedDiscreteModel{Dc}},
         fmodel::GridapDistributed.DistributedDiscreteModel{Dc},
         refinement_and_coarsening_flags::MPIArray{<:AbstractVector},
         stride) where Dc

  function setup_communication_buffers_fine_partition(cparts,
                                              fmodel,
                                              cmodel::Union{Nothing,GridapDistributed.DistributedDiscreteModel},
                                              fine_to_coarse_faces_map::Union{Nothing,MPIArray},
                                              fcell_to_child_id::Union{Nothing,MPIArray})
    fgids=get_cell_gids(fmodel)
    map(fmodel.models,fine_to_coarse_faces_map,fcell_to_child_id) do fmodel, fine_to_coarse_faces_map, fcell_to_child_id
      if GridapDistributed.i_am_in(cparts)
        cgids        = get_cell_gids(cmodel)
        cpartition   = PArrays.getany(partition(cgids))
        _setup_communication_buffers_fine_partition(fine_to_coarse_faces_map,
                                                    fcell_to_child_id,
                                                    partition(fgids),
                                                    cpartition)
      else
        # cmodel might be distributed among less processes than fmodel
        Int[], Int[], Int[], Int[]
      end 
    end |> tuple_of_arrays
  end 
   
  function _setup_communication_buffers_fine_partition(
        fine_to_coarse_faces_map::AbstractVector{<:AbstractVector{<:Integer}},
        fcell_to_child_id::AbstractVector{<:Integer},
        fpartition::MPIArray,
        cpartition::AbstractLocalIndices)
    # Note: Reversing snd and rcv
    lids_rcv,lids_snd = map(PArrays.getany,assembly_local_indices(fpartition))
    @debug "$(MPI.Comm_rank(MPI.COMM_WORLD)): lids_rcv: $(lids_rcv)"
    @debug "$(MPI.Comm_rank(MPI.COMM_WORLD)): lids_snd: $(lids_snd)"
    cgids_data = local_to_global(cpartition)[fine_to_coarse_faces_map[end][lids_snd.data]]
    cgids_snd  = PArrays.JaggedArray(cgids_data,lids_snd.ptrs)
    cgids_rcv  = PArrays.JaggedArray(Vector{Int}(undef,length(lids_rcv.data)),lids_rcv.ptrs)
    cell_child_id_data = fcell_to_child_id[lids_snd.data]
    cell_child_id_snd  = PArrays.JaggedArray(cell_child_id_data,lids_snd.ptrs)
    fcell_child_id_rcv = PArrays.JaggedArray(Vector{Int}(undef,length(lids_rcv.data)),lids_rcv.ptrs)
    return cgids_snd, cgids_rcv, cell_child_id_snd, fcell_child_id_rcv
  end

  function _setup_communication_buffers_fine_partition(
    fine_to_coarse_faces_map::AbstractVector{<:Gridap.Arrays.Table},
    fcell_to_child_id::Gridap.Arrays.Table,
    fpartition::MPIArray,
    cpartition::AbstractLocalIndices)

    # Note: Reversing snd and rcv
    lids_rcv,lids_snd = map(PArrays.getany,assembly_local_indices(fpartition))

    cgids_data  = fill(-1,length(lids_snd.data))
    fcell_child_id_data = fill(-1,length(lids_snd.data))
    
    f2c_map_ptrs = fine_to_coarse_faces_map[end].ptrs 
    f2c_map_data = fine_to_coarse_faces_map[end].data
    cl2g         = local_to_global(cpartition)
    for (i,fcell) in enumerate(lids_snd.data)
      # refinement or do nothing
      if (f2c_map_ptrs[fcell+1]-f2c_map_ptrs[fcell]==1)
        ccell = f2c_map_data[f2c_map_ptrs[fcell]]
        cgids_data[i] = cl2g[ccell]
        fcell_child_id_data[i] = fcell_to_child_id.data[f2c_map_ptrs[fcell]]
      end
    end 
    cgids_snd = PArrays.JaggedArray(cgids_data,lids_snd.ptrs)
    cgids_rcv = PArrays.JaggedArray(Vector{Int}(undef,length(lids_rcv.data)),lids_rcv.ptrs)
    cell_child_id_snd = PArrays.JaggedArray(fcell_child_id_data,lids_snd.ptrs)
    cell_child_id_rcv = PArrays.JaggedArray(Vector{Int}(undef,length(lids_rcv.data)),lids_rcv.ptrs)
    return cgids_snd, cgids_rcv, cell_child_id_snd, cell_child_id_rcv
  end

  function setup_communication_buffers_coarse_partition(cparts,
                                              cmodel::Union{Nothing,GridapDistributed.DistributedDiscreteModel},
                                              fmodel,
                                              fine_to_coarse_faces_map::Union{Nothing,MPIArray},
                                              fcell_to_child_id::Union{Nothing,MPIArray})
    fgids = get_cell_gids(fmodel)
    map(fmodel.models,fine_to_coarse_faces_map,fcell_to_child_id) do fmodel, fine_to_coarse_faces_map, fcell_to_child_id
      if (!(GridapDistributed.i_am_in(cparts)))
        # cmodel might be distributed among less processes than fmodel
        Int[], Int[], Int[], Int[]
      else
        cgids = get_cell_gids(cmodel)
        _setup_communication_buffers_coarse_partition(fine_to_coarse_faces_map,
                                                      fcell_to_child_id,
                                                      partition(fgids),
                                                      partition(cgids))
      end 
    end |> tuple_of_arrays
  end 

  function _setup_communication_buffers_coarse_partition(
    fine_to_coarse_faces_map::AbstractVector{<:Gridap.Arrays.Table},
    fcell_to_child_id::Gridap.Arrays.Table,
    fpartition::MPIArray,
    cpartition::MPIArray)

    # Note: Reversing snd and rcv
    flids_rcv,flids_snd = map(PArrays.getany,assembly_local_indices(fpartition))
    clids_rcv,clids_snd = map(PArrays.getany,assembly_local_indices(cpartition))

    fgids_data   = fill(-1,length(clids_snd.data))
    child_ids_data = fill(-1,length(clids_snd.data))

    f2c_map_ptrs   = fine_to_coarse_faces_map[end].ptrs 
    f2c_map_data   = fine_to_coarse_faces_map[end].data
    fchild_id_data = fcell_to_child_id.data
    fl2g           = local_to_global(PArrays.getany(fpartition))

    for i = 1:length(flids_snd.ptrs)-1
      for j = flids_snd.ptrs[i]:flids_snd.ptrs[i+1]-1
        fcell = flids_snd.data[j]
        # fcell coarsened ...
        if (f2c_map_ptrs[fcell+1]-f2c_map_ptrs[fcell]>1)
          for k = f2c_map_ptrs[fcell]:f2c_map_ptrs[fcell+1]-1
              ccell = f2c_map_data[k]
              child_id = fchild_id_data[k]
              # find ccell in clids_snd[i]:clids_snd[i+1]-1
              for l = clids_snd.ptrs[i]:clids_snd.ptrs[i+1]-1
                if (clids_snd.data[l]==ccell)
                  fgids_data[l]  = fl2g[fcell]
                  child_ids_data[l] = child_id 
                end
              end
          end
        end   
      end 
    end
    fgids_snd    = PArrays.JaggedArray(fgids_data,clids_snd.ptrs)
    fgids_rcv    = PArrays.JaggedArray(Vector{Int}(undef,length(clids_rcv.data)),clids_rcv.ptrs)
    child_id_snd = PArrays.JaggedArray(child_ids_data,clids_snd.ptrs)
    child_id_rcv = PArrays.JaggedArray(Vector{Int}(undef,length(clids_rcv.data)),clids_rcv.ptrs)
    return fgids_snd, fgids_rcv, child_id_snd, child_id_rcv
  end

  function _check_if_coarsen(pXest_type,
                             pXest_refinement_rule_type,
                             cpartition,
                             flags,
                             stride)
    num_children = get_num_children(pXest_type,pXest_refinement_rule_type)
    num_o_c_cells = own_length(cpartition)
    cell = 1
    while cell <= num_o_c_cells
      if (flags[cell]==coarsen_flag)
        cell,coarsen = _move_fwd_and_check_if_all_children_coarsened(flags,num_o_c_cells,cell,num_children,stride)
        if coarsen
          return true
        end 
      else 
        cell = cell+1
      end
    end
    return false
  end


  function _update_fine_to_coarse_faces_map_ptrs_with_ghost_data!(
                                      fine_to_coarse_faces_map::AbstractVector{<:Vector{<:Gridap.Arrays.Table}},
                                      fcell_to_child_id::AbstractVector{<:Gridap.Arrays.Table},
                                      cparts,
                                      cmodel,
                                      fmodel)
      # IMPORTANT NOTE: in this function we are assuming that the ptrs member variable of 
      # fine_to_coarse_face_map[] is shared by fcell_to_child_id[]. Thus, we do not need 
      # to update the latter. 

      fpartition = get_cell_gids(fmodel)
      flids_rcv,flids_snd = assembly_local_indices(partition(fpartition))
      snd_buffer,rcv_buffer=map(fine_to_coarse_faces_map,
                                fcell_to_child_id,
                                flids_rcv,
                                flids_snd) do fine_to_coarse_faces_map,fcell_to_child_id,flids_rcv,flids_snd
        if (GridapDistributed.i_am_in(cparts))
          @assert fine_to_coarse_faces_map[end].ptrs === fcell_to_child_id.ptrs
          cpartition = get_cell_gids(cmodel)
          clids_rcv,clids_snd = map(PArrays.getany,assembly_local_indices(partition(cpartition)))
          f2c_map_ptrs = fine_to_coarse_faces_map[end].ptrs
          f2c_map_data = fine_to_coarse_faces_map[end].data
          snd_buffer_data = Vector{Int}(undef,length(flids_snd.data))
          snd_buffer_data .= 0 

          for i = 1:length(flids_snd.ptrs)-1
            for j = flids_snd.ptrs[i]:flids_snd.ptrs[i+1]-1
              fcell = flids_snd.data[j]
              # fcell coarsened ...
              if (f2c_map_ptrs[fcell+1]-f2c_map_ptrs[fcell]>1)
                for k = f2c_map_ptrs[fcell]:f2c_map_ptrs[fcell+1]-1
                   ccell = f2c_map_data[k]
                   # find ccell in clids_snd[i]:clids_snd[i+1]-1
                   for l = clids_snd.ptrs[i]:clids_snd.ptrs[i+1]-1
                      if (clids_snd.data[l]==ccell)
                        snd_buffer_data[j]+=1 
                        break
                      end
                   end
                end  
              else 
                snd_buffer_data[j]=1
              end  
            end 
          end
          snd_buffer = PArrays.JaggedArray(snd_buffer_data,flids_snd.ptrs)
          rcv_buffer = PArrays.JaggedArray(Vector{Int}(undef,length(flids_rcv.data)),flids_rcv.ptrs)
          snd_buffer, rcv_buffer
        else 
          Gridap.Helpers.@notimplemented  
        end
      end |> tuple_of_arrays
      parts_rcv, parts_snd = assembly_neighbors(partition(fpartition))
      PArrays.exchange_fetch!(rcv_buffer,snd_buffer,ExchangeGraph(parts_snd,parts_rcv))
      map(rcv_buffer, fine_to_coarse_faces_map, fcell_to_child_id, partition(fpartition)) do rcv_buffer, fine_to_coarse_faces_map, fcell_to_child_id, findices
        flids_rcv,_ = map(PArrays.getany,assembly_local_indices(partition(fpartition)))
        for (i,fcell) in enumerate(flids_rcv.data)
          fine_to_coarse_faces_map[end].ptrs[fcell+1]=rcv_buffer.data[i] 
        end
        for fcell=own_length(findices)+1:local_length(findices)
          fine_to_coarse_faces_map[end].ptrs[fcell+1]=
              fine_to_coarse_faces_map[end].ptrs[fcell]+fine_to_coarse_faces_map[end].ptrs[fcell+1]
        end
        resize!(fine_to_coarse_faces_map[end].data,fine_to_coarse_faces_map[end].ptrs[end]-1)
        resize!(fcell_to_child_id.data,fine_to_coarse_faces_map[end].ptrs[end]-1)
      end
  end


  function _update_fine_to_coarse_faces_map_with_cgids_ghost_data!(
                      fine_to_coarse_faces_map::AbstractVector{<:Vector{<:Vector{<:Integer}}},
                      fcell_to_child_id::AbstractVector{<:Vector{<:Integer}},
                      cparts,
                      cmodel,
                      fgids,
                      cgids_rcv,
                      child_id_rcv)
      map(fine_to_coarse_faces_map,fcell_to_child_id,cgids_rcv,child_id_rcv) do fine_to_coarse_faces_map, 
                                                                   fcell_to_child_id, 
                                                                   cgids_rcv,
                                                                   child_id_rcv
        if (GridapDistributed.i_am_in(cparts))
          cgids = get_cell_gids(cmodel)
          cpartition = PArrays.getany(partition(cgids))
          # Note: Reversing snd and rcv
          lids_rcv,_ = map(PArrays.getany,assembly_local_indices(partition(fgids)))
          glo_to_loc = global_to_local(cpartition)
          for i in 1:length(lids_rcv.data)
            lid = lids_rcv.data[i]
            gid = cgids_rcv.data[i]
            fine_to_coarse_faces_map[end][lid] = glo_to_loc[gid]
            fcell_to_child_id[lid] = child_id_rcv.data[i]
            @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))] fcell_to_child_id.data[$(lid)]=$(child_id_rcv.data[i])"
          end
        end
      end
  end 


  function _update_fine_to_coarse_faces_map_with_cgids_ghost_data!(
                                      fine_to_coarse_faces_map::AbstractVector{<:Vector{<:Gridap.Arrays.Table}},
                                      fcell_to_child_id::AbstractVector{<:Gridap.Arrays.Table},
                                      cparts,
                                      cmodel,
                                      fgids,
                                      cgids_rcv,
                                      child_id_rcv)
      map(fine_to_coarse_faces_map,fcell_to_child_id,cgids_rcv,child_id_rcv) do fine_to_coarse_faces_map, 
                                                                                fcell_to_child_id, 
                                                                                cgids_rcv,
                                                                                child_id_rcv
        if (GridapDistributed.i_am_in(cparts))
          data  = fine_to_coarse_faces_map[end].data
          ptrs  = fine_to_coarse_faces_map[end].ptrs
          cgids = get_cell_gids(cmodel)
          cpartition = PArrays.getany(partition(cgids))
          # Note: Reversing snd and rcv
          lids_rcv,_ = map(PArrays.getany,assembly_local_indices(partition(fgids)))
          glo_to_loc = global_to_local(cpartition)
          for i in 1:length(lids_rcv.data)
            lid = lids_rcv.data[i]
            gid = cgids_rcv.data[i]
            if gid != -1
               @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))] gid=$(gid) $(ptrs[lid+1]-ptrs[lid])"
               @assert (ptrs[lid+1]-ptrs[lid])==1
               data[ptrs[lid]] = glo_to_loc[gid]
               fcell_to_child_id.data[ptrs[lid]] = child_id_rcv.data[i]
               @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))] fcell_to_child_id.data[$(ptrs[lid])]=$(child_id_rcv.data[i])"
            end
          end
        end
      end
  end

  function _update_fine_to_coarse_faces_map_with_fgids_ghost_data!(
                                      fine_to_coarse_faces_map::AbstractVector{<:Vector{<:Gridap.Arrays.Table}},
                                      fcell_to_child_id::AbstractVector{<:Gridap.Arrays.Table},
                                      cparts,
                                      cmodel,
                                      fgids,
                                      fgids_rcv,
                                      child_id_rcv)
      map(fine_to_coarse_faces_map,
          fcell_to_child_id, 
          fgids_rcv, 
          child_id_rcv, 
          partition(fgids)) do fine_to_coarse_faces_map, fcell_to_child_id, fgids_rcv, child_id_rcv, fpartition
        if (GridapDistributed.i_am_in(cparts))
          f2c_data = fine_to_coarse_faces_map[end].data
          f2c_ptrs = fine_to_coarse_faces_map[end].ptrs
          fchild_id_data = fcell_to_child_id.data
          cgids = get_cell_gids(cmodel)
          cpartition = PArrays.getany(partition(cgids))
          # Note: Reversing snd and rcv
          lids_rcv,_ = map(PArrays.getany,assembly_local_indices(partition(cgids)))
          glo_to_loc = global_to_local(fpartition)
          tmp_ptrs=Dict((i=>f2c_ptrs[i] for i=own_length(fpartition)+1:local_length(fpartition)))
          @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))] child_id_rcv.data=$(child_id_rcv.data)"
          @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))] fgids_rcv_data=$(fgids_rcv.data)"
          @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))] tmp_ptrs=$(tmp_ptrs)"
          @assert length(lids_rcv.data)==length(fgids_rcv.data)
          for i in 1:length(lids_rcv.data)
            gid_fine        = fgids_rcv.data[i]
            child_id_coarse = child_id_rcv.data[i]
            if gid_fine!=-1
              @assert child_id_coarse!=-1
              lid_coarse = lids_rcv.data[i]
              lid_fine = glo_to_loc[gid_fine]
              # Add lid_fine to lid_coarse
              pos=tmp_ptrs[lid_fine]
              f2c_data[pos]=lid_coarse
              fchild_id_data[pos]=child_id_coarse
              @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))] f2c_data[$(pos)]=$(lid_coarse))"
              @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))] fchild_id_data[$(pos)]=$(child_id_coarse))"
              tmp_ptrs[lid_fine]+=1
            end
          end
        end
      end
  end

  Gridap.Helpers.@notimplementedif cmodel==nothing

  # Check if there is at least one part in which a cell was coarsened
  cgids = get_cell_gids(cmodel)
  coarsen_array = map(partition(cgids),refinement_and_coarsening_flags) do cpartition,flags
    _check_if_coarsen(pXest_type,pXest_refinement_rule_type,cpartition,flags,stride) 
  end
  or_func(a,b)=a || b
  coarsen = PArrays.getany(reduction(or_func,coarsen_array;destination=:all,init=false))

  # Fill data for owned (from coarse cells owned by the processor)
  fgids = get_cell_gids(fmodel)
  f1,f2,f3 = map(fmodel.models,
                 partition(fgids),
                 refinement_and_coarsening_flags) do fmodel, fpartition, flags
    if (!(GridapDistributed.i_am_in(cparts)))
      # cmodel might be distributed among less processes than fmodel
      nothing, nothing, nothing
    else
      cgids        = get_cell_gids(cmodel)
      cmodel_local = PArrays.getany(cmodel.models)
      cpartition   = PArrays.getany(partition(cgids))
      fine_to_coarse_faces_map,
        fine_to_coarse_faces_dim,
          fcell_to_child_id =
          _process_owned_cells_fine_to_coarse_model_glue(pXest_type,
                                                         pXest_refinement_rule_type,
                                                         cmodel_local,
                                                         fmodel,
                                                         cpartition,
                                                         fpartition,
                                                         flags,
                                                         coarsen,
                                                         stride)
    end
  end |> tuple_of_arrays

  cgids_snd, cgids_rcv, fchild_id_snd, fchild_id_rcv = setup_communication_buffers_fine_partition(cparts,
                                                                    fmodel,
                                                                    cmodel,
                                                                    f1,
                                                                    f3)

  if coarsen 
    _update_fine_to_coarse_faces_map_ptrs_with_ghost_data!(f1,
                                                           f3,
                                                           cparts,
                                                           cmodel,
                                                           fmodel)                                                     

    fgids_snd, fgids_rcv, cchild_id_snd, cchild_id_rcv = setup_communication_buffers_coarse_partition(cparts,
                                                                                                      cmodel,
                                                                                                      fmodel,
                                                                                                      f1,
                                                                                                      f3)         
  end

  cgids = get_cell_gids(cmodel)
  cache = fetch_vector_ghost_values_cache(refinement_and_coarsening_flags,partition(cgids))
  fetch_vector_ghost_values!(refinement_and_coarsening_flags,cache) |> wait
  
  # Note: Reversing snd and rcv 
  fparts_rcv, fparts_snd = assembly_neighbors(partition(fgids))
  graph = ExchangeGraph(fparts_snd,fparts_rcv)
  PArrays.exchange_fetch!(cgids_rcv,cgids_snd,graph)
  PArrays.exchange_fetch!(fchild_id_rcv,fchild_id_snd,graph)

  if coarsen 
    cparts_rcv, cparts_snd = assembly_neighbors(partition(cgids))
    graph = ExchangeGraph(cparts_snd,cparts_rcv)
    PArrays.exchange_fetch!(fgids_rcv,fgids_snd,graph)
    PArrays.exchange_fetch!(cchild_id_rcv,cchild_id_snd,graph)
    _update_fine_to_coarse_faces_map_with_fgids_ghost_data!(f1,
                                                            f3,
                                                            cparts,
                                                            cmodel,
                                                            get_cell_gids(fmodel),
                                                            fgids_rcv,
                                                            cchild_id_rcv)
  end

  _update_fine_to_coarse_faces_map_with_cgids_ghost_data!(f1,
                                                          f3,
                                                          cparts,
                                                          cmodel,
                                                          fgids,
                                                          cgids_rcv,
                                                          fchild_id_rcv)


  # Create distributed glue
  glue = map(f1,f2,f3,refinement_and_coarsening_flags) do fine_to_coarse_faces_map, 
                                                          fine_to_coarse_faces_dim, 
                                                          fcell_to_child_id,
                                                          flags
    if (!(GridapDistributed.i_am_in(cparts)))
      nothing
    else
      ## The following lines are a replacement for WhiteRefinementRule()
      ## to have the types of rrule_nothing_flag and rrule_refinement_flag
      ## to be 100% equivalent for all type parameters
      polytope  = (Dc==2 ? QUAD : HEX)
      partition = Gridap.ReferenceFEs.tfill(1,Val{Dc}())
      ref_grid  = UnstructuredGrid(compute_reference_grid(polytope,partition))
      rrule_nothing_flag = 
         Gridap.Adaptivity.RefinementRule(Gridap.Adaptivity.WithoutRefinement(),polytope,ref_grid)
      rrule_refinement_flag = get_refinement_rule(pXest_type,pXest_refinement_rule_type)
      coarse_cell_to_rrule  = map(x -> (x==nothing_flag) ? 1 : 2,flags)
      rrules = Gridap.Arrays.CompressedArray([rrule_nothing_flag,rrule_refinement_flag],coarse_cell_to_rrule)
      @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))] fine_to_coarse_faces_map[end]: $(fine_to_coarse_faces_map[end])"
      @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))] fcell_to_child_id: $(fcell_to_child_id)"
      GT = isa(fine_to_coarse_faces_map,Vector{<:AbstractVector{<:Integer}}) ? Gridap.Adaptivity.RefinementGlue() : Gridap.Adaptivity.MixedGlue()
      AdaptivityGlue(GT,fine_to_coarse_faces_map,fcell_to_child_id,rrules)
    end
  end
  return glue
end

function _process_owned_cells_fine_to_coarse_model_glue(pXest_type,
                                                        pXest_refinement_rule_type,
                                                        cmodel::DiscreteModel{Dc},
                                                        fmodel::DiscreteModel{Dc},
                                                        cpartition,
                                                        fpartition,
                                                        flags,
                                                        coarsen,
                                                        stride) where Dc

  function _setup_fine_to_coarse_faces_map_table(pXest_type,
                                                 pXest_refinement_rule_type,
                                                 flags,
                                                 num_o_c_cells,
                                                 num_f_cells,
                                                 stride)
    num_children = get_num_children(pXest_type, pXest_refinement_rule_type)

    # Count cell children: 
    fine_to_coarse_faces_map_ptrs = Vector{Int}(undef,num_f_cells+1)
    fine_to_coarse_faces_map_ptrs[1] = 1
    old_cell = 1
    new_cell = 1
    while old_cell <= num_o_c_cells # For each coarse cell
      if (flags[old_cell]==nothing_flag)    # Cell not touched 
        fine_to_coarse_faces_map_ptrs[new_cell+1] = fine_to_coarse_faces_map_ptrs[new_cell]+1
        old_cell += 1
        new_cell += 1
      elseif (flags[old_cell]==refine_flag) # Cell is refined 
        for child = 1:num_children
          for j=1:stride
            fine_to_coarse_faces_map_ptrs[new_cell+1] = fine_to_coarse_faces_map_ptrs[new_cell]+1
            new_cell+=1
          end
        end
        old_cell=old_cell+stride
      else                             # Cell is coarsened
        @assert flags[old_cell]==coarsen_flag
        cell_fwd,coarsen = _move_fwd_and_check_if_all_children_coarsened(flags,num_o_c_cells,old_cell,num_children,stride)
        if coarsen
          @assert (cell_fwd-old_cell)==(num_children*stride)
            for j=1:stride
              fine_to_coarse_faces_map_ptrs[new_cell+1] = fine_to_coarse_faces_map_ptrs[new_cell]+num_children
              old_cell = old_cell+num_children
              new_cell = new_cell+1
            end
        else
          Gridap.Helpers.@unreachable
          # for j = new_cell:new_cell+(cell_fwd-old_cell+1) 
          #   fine_to_coarse_faces_map_ptrs[j+1] = fine_to_coarse_faces_map_ptrs[j]+1
          #   new_cell = new_cell+1
          # end
          # old_cell = cell_fwd
          # new_cell = new_cell+cell_fwd-old_cell
        end
      end
    end

    # Counts to pointers: 
    for j = new_cell:num_f_cells
      fine_to_coarse_faces_map_ptrs[j+1] = fine_to_coarse_faces_map_ptrs[j]
    end

    # Fill children data: 
    fine_to_coarse_faces_map_data = Vector{Int}(undef,fine_to_coarse_faces_map_ptrs[end]-1)
    fcell_to_child_id_data = Vector{Int}(undef,fine_to_coarse_faces_map_ptrs[end]-1)
    fcell_to_child_id_data .= -1
    old_cell = 1
    new_cell = 1
    while old_cell <= num_o_c_cells
      if (flags[old_cell]==refine_flag)      # Cell is refined
        for child = 1:num_children
          current_cell=old_cell
          for j=1:stride
            fine_to_coarse_faces_map_data[fine_to_coarse_faces_map_ptrs[new_cell]] = current_cell
            fcell_to_child_id_data[fine_to_coarse_faces_map_ptrs[new_cell]] = child
            new_cell+=1
            current_cell+=1
          end
        end
        old_cell=old_cell+stride
      elseif (flags[old_cell]==nothing_flag) # Cell is not touched
        fine_to_coarse_faces_map_data[fine_to_coarse_faces_map_ptrs[new_cell]]=old_cell
        fcell_to_child_id_data[fine_to_coarse_faces_map_ptrs[new_cell]]=1
        new_cell+=1
        old_cell+=1
      else                               # Cell is coarsened
        @assert flags[old_cell]==coarsen_flag
        cell_fwd,coarsen=_move_fwd_and_check_if_all_children_coarsened(flags,num_o_c_cells,old_cell,num_children,stride)
        if coarsen
          for child = 1:num_children
            current_cell=new_cell
            for j=1:stride
               fine_to_coarse_faces_map_data[fine_to_coarse_faces_map_ptrs[current_cell]+child-1] = old_cell
               fcell_to_child_id_data[fine_to_coarse_faces_map_ptrs[current_cell]+child-1] = child
               old_cell+=1
               current_cell+=1
            end
          end 
          new_cell=new_cell+stride
        else 
          Gridap.Helpers.@unreachable
          # for j = cell:cell_fwd-1 
          #   fine_to_coarse_faces_map_data[fine_to_coarse_faces_map_ptrs[c]]=j
          #   fcell_to_child_id_data[fine_to_coarse_faces_map_ptrs[c]]=1
          #   c = c+1
          #   cell = cell+1
          # end
        end
      end
    end
    Gridap.Arrays.Table(fine_to_coarse_faces_map_data, fine_to_coarse_faces_map_ptrs),
             Gridap.Arrays.Table(fcell_to_child_id_data, fine_to_coarse_faces_map_ptrs)
  end

  function _setup_fine_to_coarse_faces_map_vector!(pXest_type,
                                                   pXest_refinement_rule_type,
                                                   fine_to_coarse_faces_map,
                                                   fcell_to_child_id,
                                                   flags,
                                                   num_o_c_cells,
                                                   stride)

    @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))] num_o_c_cells: $(num_o_c_cells)"
    @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))] flags: $(flags)"
    @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))] len(fcell_to_child_id): $(length(fcell_to_child_id))"
    @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))] stride: $(stride)"


    # Go over all cells of coarse grid portion
    num_children = get_num_children(pXest_type,pXest_refinement_rule_type)
    c = 1
    cell=1
    while cell<=num_o_c_cells
      if flags[cell]==refine_flag
        for child = 1:num_children
          current_cell=cell
          for j=1:stride
            fine_to_coarse_faces_map[c] = current_cell
            fcell_to_child_id[c] = child
            @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))] fcell_to_child_id[$(c)]: $(child)"
            c+=1
            current_cell+=1
          end
        end
        cell=cell+stride
      elseif (flags[cell]==nothing_flag)
        fine_to_coarse_faces_map[c] = cell
        fcell_to_child_id[c] = 1
        @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))] fcell_to_child_id[$(c)]: 1"
        c+=1
        cell+=1
      elseif (flags[cell]==coarsen_flag)
        Gridap.Helpers.@notimplemented
      else
        @assert flags[cell]!=coarsen_flag 
        error("Unknown AMR flag")
      end
    end
  end
  
  num_f_cells   = num_cells(fmodel)       # Number of fine cells (owned+ghost)
  num_o_c_cells = own_length(cpartition)  # Number of coarse cells (owned)
  @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))] own_length(fpartition): $(own_length(fpartition))"
  ftopology = Gridap.Geometry.get_grid_topology(fmodel)
  if coarsen
    fine_to_coarse_faces_map = Vector{Gridap.Arrays.Table{Int,Vector{Int},Vector{Int}}}(undef,Dc+1)
    a,b = _setup_fine_to_coarse_faces_map_table(pXest_type,
                                                pXest_refinement_rule_type,
                                                flags,
                                                num_o_c_cells,
                                                num_f_cells,
                                                stride)
    fine_to_coarse_faces_map[Dc+1] = a
    fcell_to_child_id = b
    # In the future we should also have here the code to also setup
    # fine_to_coarse_faces_map for lower dimensional objects
  else
    fine_to_coarse_faces_map = Vector{Vector{Int}}(undef,Dc+1)
    fine_to_coarse_faces_map[Dc+1] = Vector{Int}(undef,num_f_cells)
    fcell_to_child_id = Vector{Int}(undef,num_f_cells)
    fcell_to_child_id .= -1
    _setup_fine_to_coarse_faces_map_vector!(pXest_type,
                                            pXest_refinement_rule_type,
                                            fine_to_coarse_faces_map[Dc+1],
                                            fcell_to_child_id,flags,
                                            num_o_c_cells,
                                            stride)
    for d=1:Dc
      fine_to_coarse_faces_map[d] = Vector{Int}(undef,Gridap.Geometry.num_faces(ftopology,d-1))
    end
    # Code commented out as not currently required. 
    # c_cell_faces=[]
    # cache_c_cell_faces=[]
    # f_cell_faces=[]
    # cache_f_cell_faces=[]
    # for d=1:Dc
    #   push!(c_cell_faces, Gridap.Geometry.get_faces(ctopology, Dc, d-1))
    #   push!(cache_c_cell_faces,array_cache(last(c_cell_faces)))
    #   push!(f_cell_faces, Gridap.Geometry.get_faces(ftopology, Dc, d-1))
    #   push!(cache_f_cell_faces,array_cache(last(f_cell_faces)))
    # end
    # parent_cell_faces=Vector{Vector{Int}}(undef,Dc)
    # for cell=1:num_f_cells
    #   parent_cell=fine_to_coarse_faces_map[Dc+1][cell]
    #   child=fcell_to_child_id[cell]
    #   for d=1:Dc
    #     parent_cell_faces[d]=getindex!(cache_c_cell_faces[d],
    #                                   c_cell_faces[d],
    #                                   parent_cell)
    #   end
    #   for d=1:Dc
    #     cell_f_faces=getindex!(cache_f_cell_faces[d],
    #                             f_cell_faces[d],
    #                             cell)
    #     for (lf,f) in enumerate(cell_f_faces)
    #       c     = rrule_f_to_c_lid_2D[d][child][lf]
    #       dim_c = rrule_f_to_c_dim_2D[d][child][lf]
    #       if (dim_c == Dc)
    #         fine_to_coarse_faces_map[d][f]=parent_cell
    #       else
    #         fine_to_coarse_faces_map[d][f]=parent_cell_faces[dim_c+1][c]
    #       end
    #       fine_to_coarse_faces_dim[d][f]=dim_c
    #     end
    #   end
    # end
  end
  # To-think: how fine_to_coarse_faces_dim should look like whenever we have 
  # coarsening and refinement at the same time? By now it is not being used,
  # so least concern. 
  fine_to_coarse_faces_dim = Vector{Vector{Int}}(undef,Dc)
  for d=1:Dc
    fine_to_coarse_faces_dim[d] = Vector{Int}(undef,Gridap.Geometry.num_faces(ftopology,d-1))
  end
  return fine_to_coarse_faces_map, fine_to_coarse_faces_dim, fcell_to_child_id
end


function Gridap.Adaptivity.refine(model::OctreeDistributedDiscreteModel{Dc,Dp}; parts=nothing) where {Dc,Dp}
   old_comm = model.parts.comm
   if (GridapDistributed.i_am_in(old_comm))
     # Copy and refine input p4est
     ptr_new_pXest = pXest_copy(model.pXest_type, model.ptr_pXest)
     pXest_uniformly_refine!(model.pXest_type, ptr_new_pXest)
   else
     ptr_new_pXest = nothing
   end

   new_comm = isa(parts,Nothing) ? old_comm : parts.comm
   if GridapDistributed.i_am_in(new_comm)
      pXest_type = _dim_to_pXest_type(Dc)

      if !isa(parts,Nothing)
        aux = ptr_new_pXest
        ptr_new_pXest = _pXest_to_new_comm(pXest_type,
                                           ptr_new_pXest,
                                           model.ptr_pXest_connectivity,
                                           model.parts.comm,
                                           parts.comm)
        if GridapDistributed.i_am_in(old_comm)
          pXest_destroy(pXest_type,aux)
        end
      end

      # Extract ghost and lnodes
      ptr_pXest_ghost  = setup_pXest_ghost(pXest_type, ptr_new_pXest)
      ptr_pXest_lnodes = setup_pXest_lnodes(pXest_type, ptr_new_pXest, ptr_pXest_ghost)

      # Build fine-grid mesh
      new_parts = isa(parts,Nothing) ? model.parts : parts
      fmodel = setup_distributed_discrete_model(pXest_type,
                                              new_parts,
                                              model.coarse_model,
                                              model.ptr_pXest_connectivity,
                                              ptr_new_pXest,
                                              ptr_pXest_ghost,
                                              ptr_pXest_lnodes)

      pXest_lnodes_destroy(pXest_type,ptr_pXest_lnodes)
      pXest_ghost_destroy(pXest_type,ptr_pXest_ghost)

      dglue = _compute_fine_to_coarse_model_glue(model.parts,
                                                 model.dmodel,
                                                 fmodel)

      non_conforming_glue = _create_conforming_model_non_conforming_glue(fmodel)
                                           

      ref_model = OctreeDistributedDiscreteModel(Dc,Dp,
                                     new_parts,
                                     fmodel,
                                     non_conforming_glue,
                                     model.coarse_model,
                                     model.ptr_pXest_connectivity,
                                     ptr_new_pXest,
                                     pXest_type,
                                     PXestUniformRefinementRuleType(),
                                     false,
                                     model)

      return ref_model, dglue
   else
    new_parts = isa(parts,Nothing) ? model.parts : parts
    return VoidOctreeDistributedDiscreteModel(model,new_parts), nothing
   end
end


function _refine_coarsen_balance!(model::OctreeDistributedDiscreteModel{Dc,Dp}, 
                                  refinement_and_coarsening_flags::MPIArray{<:Vector}) where {Dc,Dp}

  
  pXest_type = model.pXest_type
  init_fn_callback_c = pXest_reset_callbacks(pXest_type)
  coarsen_fn_callback_c = pXest_coarsen_callbacks(pXest_type)
  refine_callback_c,refine_replace_callback_c = pXest_refine_callbacks(pXest_type)

  map(model.dmodel.models,refinement_and_coarsening_flags) do lmodel, flags
    # The length of the local flags array has to match the number of 
    # cells in the model. This includes both owned and ghost cells. 
    # Only the flags for owned cells are actually taken into account. 
    @assert num_cells(lmodel)==length(flags)
    pXest_reset_data!(pXest_type, model.ptr_pXest, Cint(sizeof(Cint)), init_fn_callback_c, pointer(flags))
  end

  # Copy input p4est, refine and balance
  ptr_new_pXest = pXest_copy(pXest_type, model.ptr_pXest)
  pXest_refine!(pXest_type, ptr_new_pXest,
                refine_callback_c,
                refine_replace_callback_c)
  pXest_coarsen!(pXest_type, ptr_new_pXest, coarsen_fn_callback_c)
  pXest_balance!(pXest_type, ptr_new_pXest)
  pXest_update_flags!(pXest_type,model.ptr_pXest,ptr_new_pXest)
  ptr_new_pXest
end 

function Gridap.Adaptivity.adapt(model::OctreeDistributedDiscreteModel{Dc,Dp}, 
		                             refinement_and_coarsening_flags::MPIArray{<:Vector{<:Integer}};
                                 parts=nothing) where {Dc,Dp}

  Gridap.Helpers.@notimplementedif parts!=nothing

  _refinement_and_coarsening_flags = map(refinement_and_coarsening_flags) do flags
    convert(Vector{Cint},flags)
   end 
  
  ptr_new_pXest = _refine_coarsen_balance!(model, _refinement_and_coarsening_flags)

  # Extract ghost and lnodes
  ptr_pXest_ghost  = setup_pXest_ghost(model.pXest_type, ptr_new_pXest)
  ptr_pXest_lnodes = setup_pXest_lnodes_nonconforming(model.pXest_type, ptr_new_pXest, ptr_pXest_ghost)

  # Build fine-grid mesh
  fmodel,non_conforming_glue = setup_non_conforming_distributed_discrete_model(model.pXest_type,
                                                                              model.parts,
                                                                              model.coarse_model,
                                                                              model.ptr_pXest_connectivity,
                                                                              ptr_new_pXest,
                                                                              ptr_pXest_ghost,
                                                                              ptr_pXest_lnodes)
    
  pXest_ghost_destroy(model.pXest_type,ptr_pXest_ghost)
  pXest_lnodes_destroy(model.pXest_type,ptr_pXest_lnodes)
  pXest_refinement_rule_type = PXestUniformRefinementRuleType()
  stride = pXest_stride_among_children(model.pXest_type,pXest_refinement_rule_type,model.ptr_pXest)
  adaptivity_glue = _compute_fine_to_coarse_model_glue(model.pXest_type,
                                                       pXest_refinement_rule_type,
                                                       model.parts,
                                                       model.dmodel,
                                                       fmodel,
                                                       _refinement_and_coarsening_flags,
                                                       stride)
  adaptive_models = map(local_views(model),
                        local_views(fmodel),
                        adaptivity_glue) do model, fmodel, glue 
      Gridap.Adaptivity.AdaptedDiscreteModel(fmodel,model,glue)
  end
  fmodel = GridapDistributed.GenericDistributedDiscreteModel(adaptive_models,get_cell_gids(fmodel))
  ref_model = OctreeDistributedDiscreteModel(Dc,Dp,
                                             model.parts,
                                             fmodel,
                                             non_conforming_glue,
                                             model.coarse_model,
                                             model.ptr_pXest_connectivity,
                                             ptr_new_pXest,
                                             model.pXest_type,
                                             pXest_refinement_rule_type,
                                             false,
                                             model)
  return ref_model, adaptivity_glue
end

function Gridap.Adaptivity.coarsen(model::OctreeDistributedDiscreteModel{Dc,Dp}) where {Dc,Dp}
  comm = model.parts.comm
  if (GridapDistributed.i_am_in(comm))
    # Copy and coarsen input p4est
    ptr_new_pXest = pXest_copy(model.pXest_type, model.ptr_pXest)
    pXest_uniformly_coarsen!(model.pXest_type, ptr_new_pXest)
  else
    ptr_new_pXest=nothing
  end

  if (GridapDistributed.i_am_in(comm))
     # Extract ghost and lnodes
     ptr_pXest_ghost  = setup_pXest_ghost(model.pXest_type, ptr_new_pXest)
     ptr_pXest_lnodes = setup_pXest_lnodes(model.pXest_type, ptr_new_pXest, ptr_pXest_ghost)

     # Build coarse-grid mesh
     cmodel = setup_distributed_discrete_model(model.pXest_type,
                                             model.parts,
                                             model.coarse_model,
                                             model.ptr_pXest_connectivity,
                                             ptr_new_pXest,
                                             ptr_pXest_ghost,
                                             ptr_pXest_lnodes)

     pXest_lnodes_destroy(model.pXest_type,ptr_pXest_lnodes)
     pXest_ghost_destroy(model.pXest_type,ptr_pXest_ghost)

     dglue = _compute_fine_to_coarse_model_glue(model.parts,
                                                cmodel,
                                                model.dmodel)

     nc_glue=_create_conforming_model_non_conforming_glue(cmodel)

     pXest_type = _dim_to_pXest_type(Dc)

     c_octree_model = OctreeDistributedDiscreteModel(Dc,Dp,
                                    model.parts,
                                    cmodel,
                                    nc_glue,
                                    model.coarse_model,
                                    model.ptr_pXest_connectivity,
                                    ptr_new_pXest,
                                    pXest_type,
                                    PXestUniformRefinementRuleType(),
                                    false,
                                    model)
     return c_octree_model, dglue
  else
     return VoidOctreeDistributedDiscreteModel(model,model.parts), nothing
  end
end



# We have a p4est distributed among P processors. This function
# instantiates the same among Q processors.
function _pXest_to_new_comm(pXest_type,ptr_pXest, ptr_pXest_conn, old_comm, new_comm)
  A = is_included(old_comm,new_comm) # old \subset new (smaller to larger nparts)
  B = is_included(new_comm,old_comm) # old \supset new (larger to smaller nparts)
  @assert xor(A,B)
  if (A)
    _pXest_to_new_comm_old_subset_new(pXest_type,ptr_pXest, ptr_pXest_conn, old_comm, new_comm)
  else
    _pXest_to_new_comm_old_supset_new(pXest_type,ptr_pXest, ptr_pXest_conn, old_comm, new_comm)
  end
end

function _pXest_to_new_comm_old_subset_new(pXest_type,ptr_pXest, ptr_pXest_conn, old_comm, new_comm)
  if (GridapDistributed.i_am_in(new_comm))
    new_comm_num_parts    = GridapDistributed.num_parts(new_comm)
    global_first_quadrant = Vector{P4est_wrapper.p4est_gloidx_t}(undef,new_comm_num_parts+1)

    pXest_conn = ptr_pXest_conn[]
    pertree = Vector{P4est_wrapper.p4est_gloidx_t}(undef,pXest_conn.num_trees+1)
    if (GridapDistributed.i_am_in(old_comm))
      pXest = ptr_pXest[]
      old_comm_num_parts = GridapDistributed.num_parts(old_comm)
      old_global_first_quadrant = unsafe_wrap(Array,
                                              pXest.global_first_quadrant,
                                              old_comm_num_parts+1)

      for i = 1:length(old_global_first_quadrant)
        global_first_quadrant[i] = old_global_first_quadrant[i]
      end
      for i = length(old_global_first_quadrant)+1:length(global_first_quadrant)
        global_first_quadrant[i] = old_global_first_quadrant[end]
      end
      MPI.Bcast!(global_first_quadrant,0,new_comm)
      quadrants = pXest_deflate_quadrants(pXest_type,ptr_pXest,C_NULL)
      pXest_comm_count_pertree(pXest_type,ptr_pXest,pertree)
      MPI.Bcast!(pertree,0,new_comm)
    else
      MPI.Bcast!(global_first_quadrant,0,new_comm)
      quadrants = sc_array_new_count(sizeof(p4est_quadrant_t), 0)
      MPI.Bcast!(pertree,0,new_comm)
    end
    return pXest_inflate(pXest_type,
                         new_comm,
                         ptr_pXest_conn,
                         global_first_quadrant,
                         pertree,
                         quadrants,
                         C_NULL,
                         C_NULL)
  else
    return nothing
  end
end

function _pXest_to_new_comm_old_supset_new(pXest_type,ptr_pXest, ptr_pXest_conn, old_comm, new_comm)
  @assert GridapDistributed.i_am_in(old_comm)
  pXest = ptr_pXest[]
  pXest_conn = ptr_pXest_conn[]

  pertree = Vector{P4est_wrapper.p4est_gloidx_t}(undef,pXest_conn.num_trees+1)
  pXest_comm_count_pertree(pXest_type,ptr_pXest,pertree)

  if (GridapDistributed.i_am_in(new_comm))
    new_comm_num_parts = GridapDistributed.num_parts(new_comm)
    global_first_quadrant = Vector{P4est_wrapper.p4est_gloidx_t}(undef,new_comm_num_parts+1)
    old_comm_num_parts = GridapDistributed.num_parts(old_comm)
    old_global_first_quadrant = unsafe_wrap(Array,
                                            pXest.global_first_quadrant,
                                            old_comm_num_parts+1)

    new_global_first_quadrant = unsafe_wrap(Array,
                                            pXest.global_first_quadrant,
                                            new_comm_num_parts+1)

    for i = 1:length(new_global_first_quadrant)
      global_first_quadrant[i] = old_global_first_quadrant[i]
    end
    quadrants = pXest_deflate_quadrants(pXest_type,ptr_pXest,C_NULL)

    return pXest_inflate(pXest_type,
                         new_comm,
                         ptr_pXest_conn,
                         global_first_quadrant,
                         pertree,
                         quadrants,
                         C_NULL,
                         C_NULL)
  else
    return nothing
  end
end


function is_included(partsA,partsB)
  @assert GridapDistributed.i_am_in(partsA.comm) || GridapDistributed.i_am_in(partsB.comm)
  is_included(partsA.comm,partsB.comm)
end

function is_included(commA::MPI.Comm,commB::MPI.Comm)
  @assert GridapDistributed.i_am_in(commA) || GridapDistributed.i_am_in(commB)
  num_partsA=GridapDistributed.num_parts(commA)
  num_partsB=GridapDistributed.num_parts(commB)
  if (num_partsA==num_partsB)
    return false
  end
  if (GridapDistributed.i_am_in(commA) && GridapDistributed.i_am_in(commB))
    result=num_partsA < num_partsB
    if (result)
      result=MPI.Allreduce(Int8(result),MPI.LOR,commB)
    else
      result=MPI.Allreduce(Int8(result),MPI.LOR,commA)
    end
  else
    result=false
    if (GridapDistributed.i_am_in(commB))
      result=MPI.Allreduce(Int8(result),MPI.LOR,commB)
    else
      result=MPI.Allreduce(Int8(result),MPI.LOR,commA)
    end
  end
  Bool(result)
end

# Assumptions. Either:
# A) model.parts MPI tasks are included in parts_redistributed_model MPI tasks; or
# B) model.parts MPI tasks include parts_redistributed_model MPI tasks
const WeightsArrayType=Union{Nothing,MPIArray{<:Vector{<:Integer}}}
function GridapDistributed.redistribute(model::OctreeDistributedDiscreteModel{Dc,Dp}, 
                                        parts_redistributed_model=model.parts;
                                        weights::WeightsArrayType=nothing) where {Dc,Dp}
  parts = (parts_redistributed_model === model.parts) ? model.parts : parts_redistributed_model
  _weights=nothing
  if (weights !== nothing)
    Gridap.Helpers.@notimplementedif parts!==model.parts
    _weights=map(model.dmodel.models,weights) do lmodel,weights
       # The length of the local weights array has to match the number of 
       # cells in the model. This includes both owned and ghost cells. 
       # Only the flags for owned cells are actually taken into account. 
       @assert num_cells(lmodel)==length(weights)
       convert(Vector{Cint},weights)
    end
  end

  comm  = parts.comm
  if (GridapDistributed.i_am_in(model.parts.comm) || GridapDistributed.i_am_in(parts.comm))
    if (parts_redistributed_model !== model.parts)
      A=is_included(model.parts,parts_redistributed_model)
      B=is_included(parts_redistributed_model,model.parts)
      @assert A || B
    end
    if (parts_redistributed_model===model.parts || A)
      _redistribute_parts_subseteq_parts_redistributed(model,parts_redistributed_model,_weights)
    else
      _redistribute_parts_supset_parts_redistributed(model, parts_redistributed_model)
    end
  else
    VoidOctreeDistributedDiscreteModel(model,model.parts), nothing
  end
end

function _redistribute_parts_subseteq_parts_redistributed(model::OctreeDistributedDiscreteModel{Dc,Dp},
                                                          parts_redistributed_model, 
                                                          _weights::WeightsArrayType) where {Dc,Dp}
  parts = (parts_redistributed_model === model.parts) ? model.parts : parts_redistributed_model
  if (parts_redistributed_model === model.parts)
    ptr_pXest_old = model.ptr_pXest
  else
    ptr_pXest_old = _pXest_to_new_comm(model.pXest_type,
                                       model.ptr_pXest,
                                       model.ptr_pXest_connectivity,
                                       model.parts.comm,
                                       parts.comm)
  end
  ptr_pXest_new = pXest_copy(model.pXest_type, ptr_pXest_old)
  if (_weights !== nothing)
    init_fn_callback_c = pXest_reset_callbacks(model.pXest_type)
    map(_weights) do _weights 
      pXest_reset_data!(model.pXest_type, ptr_pXest_new, Cint(sizeof(Cint)), init_fn_callback_c, pointer(_weights))
    end
    pXest_partition!(model.pXest_type, ptr_pXest_new; weights_set=true)
  else
    pXest_partition!(model.pXest_type, ptr_pXest_new; weights_set=false)
  end 

  # Compute RedistributeGlue
  parts_snd, lids_snd, old2new = pXest_compute_migration_control_data(model.pXest_type,ptr_pXest_old,ptr_pXest_new)
  parts_rcv, lids_rcv, new2old = pXest_compute_migration_control_data(model.pXest_type,ptr_pXest_new,ptr_pXest_old)

  lids_rcv, parts_rcv, lids_snd, parts_snd, old2new, new2old =
       _to_pdata(parts, lids_rcv, parts_rcv, lids_snd, parts_snd, old2new, new2old)

  glue = GridapDistributed.RedistributeGlue(parts,model.parts,parts_rcv,parts_snd,lids_rcv,lids_snd,old2new,new2old)

  # Extract ghost and lnodes
  ptr_pXest_ghost  = setup_pXest_ghost(model.pXest_type, ptr_pXest_new)
  ptr_pXest_lnodes = setup_pXest_lnodes_nonconforming(model.pXest_type, ptr_pXest_new, ptr_pXest_ghost)

  # Build fine-grid mesh
  fmodel, non_conforming_glue = setup_non_conforming_distributed_discrete_model(model.pXest_type,
                                            parts,
                                            model.coarse_model,
                                            model.ptr_pXest_connectivity,
                                            ptr_pXest_new,
                                            ptr_pXest_ghost,
                                            ptr_pXest_lnodes)

  pXest_lnodes_destroy(model.pXest_type,ptr_pXest_lnodes)
  pXest_ghost_destroy(model.pXest_type,ptr_pXest_ghost)

  red_model = OctreeDistributedDiscreteModel(Dc,Dp,
                                             parts,
                                             fmodel,
                                             non_conforming_glue,
                                             model.coarse_model,
                                             model.ptr_pXest_connectivity,
                                             ptr_pXest_new,
                                             model.pXest_type,
                                             model.pXest_refinement_rule_type,
                                             false,
                                             model)
  return red_model, glue
end

function _redistribute_parts_supset_parts_redistributed(
     model::OctreeDistributedDiscreteModel{Dc,Dp}, 
     parts_redistributed_model) where {Dc,Dp}
  @assert model.parts !== parts_redistributed_model

  subset_comm = parts_redistributed_model.comm
  supset_comm = model.parts.comm
  N=num_cells(model)
  if (GridapDistributed.i_am_in(subset_comm))
    # This piece of code replicates the logic behind the
    # "p4est_partition_cut_gloidx" function in the p4est library
    psub=PArrays.getany(parts_redistributed_model)
    Psub=GridapDistributed.num_parts(subset_comm)
    first_global_quadrant=Int64(floor((Float64(N)*Float64(psub-1))/(Float64(Psub))))
    @assert first_global_quadrant>=0 && first_global_quadrant<N
  else
    first_global_quadrant=N
  end

  Psup=GridapDistributed.num_parts(supset_comm)
  psub=PArrays.getany(model.parts)
  num_cells_per_part=Vector{Cint}(undef, Psup+1)
  parts_offsets=MPI.Allgather(first_global_quadrant,supset_comm)
  num_cells_per_part[1:length(parts_offsets)] .= parts_offsets
  num_cells_per_part[end]=N
  for i=1:length(num_cells_per_part)-1
    num_cells_per_part[i]=num_cells_per_part[i+1]-num_cells_per_part[i]
  end
  # p4est_vtk_write_file(model.ptr_pXest, C_NULL, "model.ptr_pXest")

  ptr_pXest_old=pXest_copy(model.pXest_type,model.ptr_pXest)
  pXest_partition_given!(model.pXest_type, ptr_pXest_old, num_cells_per_part)

  # p4est_vtk_write_file(ptr_pXest_old, C_NULL, "ptr_pXest_old")

  # ptr_pXest_old is distributed over supset_comm
  # once created, ptr_pXest_new is distributed over subset_comm
  ptr_pXest_new = _pXest_to_new_comm(model.pXest_type,
                                     ptr_pXest_old,
                                     model.ptr_pXest_connectivity,
                                     supset_comm,
                                     subset_comm)

  # Compute RedistributeGlue
  parts_snd, lids_snd, old2new =
      pXest_compute_migration_control_data(model.pXest_type,model.ptr_pXest,ptr_pXest_old)
  parts_rcv, lids_rcv, new2old =
      pXest_compute_migration_control_data(model.pXest_type,ptr_pXest_old,model.ptr_pXest)

  pXest_destroy(model.pXest_type,ptr_pXest_old)

  lids_rcv, parts_rcv, lids_snd, parts_snd, old2new, new2old =
       _to_pdata(model.parts, lids_rcv, parts_rcv, lids_snd, parts_snd, old2new, new2old)

  glue = GridapDistributed.RedistributeGlue(parts_redistributed_model,model.parts,
                                            parts_rcv,parts_snd,lids_rcv,lids_snd,
                                            old2new,new2old)

  if (GridapDistributed.i_am_in(subset_comm))
    # p4est_vtk_write_file(ptr_pXest_new, C_NULL, "ptr_pXest_new")

    pXest_type = _dim_to_pXest_type(Dc)

    # Extract ghost and lnodes
    ptr_pXest_ghost  = setup_pXest_ghost(pXest_type, ptr_pXest_new)
    ptr_pXest_lnodes = setup_pXest_lnodes(pXest_type, ptr_pXest_new, ptr_pXest_ghost)

    # # Build fine-grid mesh
    fmodel = setup_distributed_discrete_model(pXest_type,
                                              parts_redistributed_model,
                                              model.coarse_model,
                                              model.ptr_pXest_connectivity,
                                              ptr_pXest_new,
                                              ptr_pXest_ghost,
                                              ptr_pXest_lnodes)


    pXest_lnodes_destroy(pXest_type,ptr_pXest_lnodes)
    pXest_ghost_destroy(pXest_type,ptr_pXest_ghost)

    non_conforming_glue = _create_conforming_model_non_conforming_glue(fmodel)

    red_model = OctreeDistributedDiscreteModel(Dc,Dp,
                                               parts_redistributed_model,
                                               fmodel,
                                               non_conforming_glue,
                                               model.coarse_model,
                                               model.ptr_pXest_connectivity,
                                               ptr_pXest_new,
                                               pXest_type,
                                               model.pXest_refinement_rule_type,
                                               false,
                                               model)
    return red_model, glue
  else
    return VoidOctreeDistributedDiscreteModel(model,parts_redistributed_model), nothing
  end
end

function _to_pdata(parts, lids_rcv, parts_rcv, lids_snd, parts_snd, old2new, new2old)
  lids_rcv, parts_rcv = map(parts) do _
    lids_rcv, parts_rcv
  end |> tuple_of_arrays
  lids_snd, parts_snd = map(parts) do _
    lids_snd, parts_snd
  end |> tuple_of_arrays
  old2new, new2old = map(parts) do _
    old2new,new2old
  end |> tuple_of_arrays
  lids_rcv, parts_rcv, lids_snd, parts_snd, old2new, new2old
end

# In the local scope of this function, the term "face"
# should be understood as a generic d-face, i.e., 
# either a vertex, edge, face, etc. 
function process_current_face!(gridap_cell_faces,
  regular_face_p4est_to_gridap,
  num_regular_faces,
  p4est_faces,
  p4est_lface,
  p4est_gface,
  p4est_lface_to_gridap_lface)

  if !(haskey(regular_face_p4est_to_gridap, p4est_gface))
    num_regular_faces += 1
    regular_face_p4est_to_gridap[p4est_gface] = num_regular_faces
  end
  gridap_cell_faces[p4est_lface_to_gridap_lface[p4est_lface]] =
    regular_face_p4est_to_gridap[p4est_gface]
  return num_regular_faces
end

const p4est_face_corners = [0 2; 1 3; 0 1; 2 3]

const p8est_face_corners = [ 0 2 4 6 ;
                             1 3 5 7 ;
                             0 1 4 5 ;
                             2 3 6 7 ;
                             0 1 2 3 ;
                             4 5 6 7 ]
  
const p6est_subface_to_hanging_edges_within_subface = 
[ 
  1;
  0;
  3;
  2;
] 

const p8est_subface_to_hanging_edges_within_subface = 
[ 
  1 3;
  1 2;
  0 3;
  0 2;
]   

# Here, we assign the identifier 5 to a "horizontal" hanging edge
# within a face, and the identifier 6 to a "vertical" hanging edge
# within a face
const p6est_subface_to_hanging_edges_within_face = 
[ 
  5;
  5;
  6;
  6;
] 

# Here, we are assigning arbitrarily a local identifier
# within the master face to each hanging edge. 
# 1: south, 2: north, 3: west and 4: east.
const p8est_subface_to_hanging_edges_within_face = 
[ 
  3 1;
  4 1;
  3 2;
  4 2;
]   

const p8est_edge_corners = [ 0  1;
                             2  3;
                             4  5;
                             6  7;
                             0  2;
                             1  3;
                             4  6;
                             5  7;
                             0  4;
                             1  5;
                             2  6;
                             3  7 ]


const hanging_vertex_code         = -2
const hanging_edge_from_face_code = -3
function num_cell_vertices(::Type{Val{Dc}}) where Dc
  2^Dc
end 

function num_cell_edges(::Type{Val{Dc}}) where Dc 
  Dc==2 ? 0 : 12
end 

function num_cell_faces(::Type{Val{Dc}}) where Dc
  2*Dc
end 

function setup_non_conforming_distributed_discrete_model(pXest_type::PXestType,
                                                         parts,
                                                         coarse_discrete_model,
                                                         ptr_pXest_connectivity,
                                                         ptr_pXest,
                                                         ptr_pXest_ghost,
                                                         ptr_pXest_lnodes)

  Dc=num_cell_dims(pXest_type)                                                       

  cell_prange = setup_cell_prange(pXest_type, parts, ptr_pXest, ptr_pXest_ghost)

  gridap_cell_faces,
  non_conforming_glue=
    generate_cell_faces_and_non_conforming_glue(pXest_type,ptr_pXest_lnodes, cell_prange)

  cell_corner_lids = map(gridap_cell_faces[1]) do cell_lids
    Gridap.Arrays.Table(cell_lids.data,cell_lids.ptrs) # JaggedArray -> Table
  end

  cell_vertex_coordinates = generate_cell_vertex_coordinates(
    pXest_type,
    cell_corner_lids,
    ptr_pXest_connectivity,
    ptr_pXest,
    ptr_pXest_ghost
  )
  grid, topology = generate_grid_and_topology(
    pXest_type,cell_corner_lids,cell_vertex_coordinates
  )

  map(topology,gridap_cell_faces[Dc]) do topology,cell_faces
    cell_faces_gridap = Gridap.Arrays.Table(cell_faces.data,cell_faces.ptrs)
    topology.n_m_to_nface_to_mfaces[Dc+1,Dc] = cell_faces_gridap
    topology.n_m_to_nface_to_mfaces[Dc,Dc+1] = Gridap.Geometry.generate_cells_around(cell_faces_gridap)
  end

  if (Dc==3)
    map(topology,gridap_cell_faces[Dc-1]) do topology,cell_edges
      cell_edges_gridap = Gridap.Arrays.Table(cell_edges.data,cell_edges.ptrs)
      topology.n_m_to_nface_to_mfaces[Dc+1,Dc-1] = cell_edges_gridap
      topology.n_m_to_nface_to_mfaces[Dc-1,Dc+1] = Gridap.Geometry.generate_cells_around(cell_edges_gridap)
    end
  end 

  face_labeling=generate_face_labeling(pXest_type,
                                       parts,
                                       cell_prange,
                                       coarse_discrete_model,
                                       topology,
                                       ptr_pXest,
                                       ptr_pXest_ghost)

  coarse_face_labeling = get_face_labeling(coarse_discrete_model)

  _set_hanging_labels!(face_labeling,non_conforming_glue,coarse_face_labeling)

  discretemodel=map(grid,topology,face_labeling) do grid, topology, face_labeling
    Gridap.Geometry.UnstructuredDiscreteModel(grid,topology,face_labeling)
  end
  GridapDistributed.DistributedDiscreteModel(discretemodel,cell_prange), non_conforming_glue
end

function _compute_max_entity_id(face_labeling)
  max_entity_id = typemin(eltype(first(face_labeling.d_to_dface_to_entity))) 
  for i=1:length(face_labeling.d_to_dface_to_entity)
    max_entity_id=max(maximum(face_labeling.d_to_dface_to_entity[i]),max_entity_id)
  end
  max_entity_id
end 

function _compute_min_entity_id(face_labeling)
  min_entity_id = typemin(eltype(first(face_labeling.d_to_dface_to_entity))) 
  for i=1:length(face_labeling.d_to_dface_to_entity)
    min_entity_id=min(minimum(face_labeling.d_to_dface_to_entity[i]),min_entity_id)
  end
  min_entity_id
end 

function _set_hanging_labels!(face_labeling,non_conforming_glue,coarse_face_labeling)
  max_entity_id = _compute_max_entity_id(coarse_face_labeling)
  hanging_entitity_ids = Dict{Int,Bool}()
  map(face_labeling,
      non_conforming_glue) do face_labeling,ncglue 
    for i=1:length(ncglue.num_hanging_faces)
      num_regular_faces_i = ncglue.num_regular_faces[i]
      num_hanging_faces_i = ncglue.num_hanging_faces[i]
      for j=num_regular_faces_i+1:num_regular_faces_i+num_hanging_faces_i
        hanging_entity_id = max_entity_id + face_labeling.d_to_dface_to_entity[i][j]
        @debug "hanging $(i-1)-face: $(j) hanging_entity_id=$(hanging_entity_id) face_labeling.d_to_dface_to_entity[$(i)][$(j)]=$(face_labeling.d_to_dface_to_entity[i][j])"
        face_labeling.d_to_dface_to_entity[i][j]=hanging_entity_id
        hanging_entitity_ids[hanging_entity_id]=true
      end      
    end
  end
  map(face_labeling) do face_labeling 
    add_tag!(face_labeling,"hanging",collect(keys(hanging_entitity_ids)))
  end
end 




