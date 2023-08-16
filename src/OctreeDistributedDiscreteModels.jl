
const nothing_flag  = Cint(0)
const refine_flag   = Cint(1)
const coarsen_flag  = Cint(2)

struct NonConformingGlue{Dc,
                         A,
                         B,
                         C}
                         #A<:AbstractVector{<:AbstractVector{<:Integer}},
                         #B<:AbstractVector{<:AbstractVector{<:Integer}},
                         #C<:AbstractVector{<:AbstractVector{<:Tuple{<:Integer,:Integer,Integer}}}}
  num_regular_faces  :: A 
  num_hanging_faces  :: B
  hanging_faces_glue :: C
  function NonConformingGlue(num_regular_faces,
                             num_hanging_faces,
                             hanging_faces_glue)
    Dc=length(num_regular_faces)
    @assert length(num_hanging_faces)==Dc
    @assert length(hanging_faces_glue)==Dc
    A=typeof(num_regular_faces)
    B=typeof(num_hanging_faces)
    C=typeof(hanging_faces_glue)
    new{Dc,A,B,C}(num_regular_faces,
                  num_hanging_faces,
                  hanging_faces_glue)
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
    NonConformingGlue(num_regular_faces,
                      num_hanging_faces,
                      hanging_faces_glue)
  end
end

mutable struct OctreeDistributedDiscreteModel{Dc,Dp,A,B,C,D,E,F} <: GridapDistributed.DistributedDiscreteModel{Dc,Dp}
  parts                       :: A
  dmodel                      :: B
  non_conforming_glue         :: C 
  coarse_model                :: D
  ptr_pXest_connectivity      :: E
  ptr_pXest                   :: F

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
    owns_ptr_pXest_connectivity::Bool,
    gc_ref)

    if (isa(dmodel,GridapDistributed.DistributedDiscreteModel))
      Gridap.Helpers.@check Dc == Gridap.Geometry.num_cell_dims(dmodel)
      Gridap.Helpers.@check Dc == Gridap.Geometry.num_point_dims(dmodel)
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
                                        owns_ptr_pXest_connectivity,
                                        gc_ref)
end


function OctreeDistributedDiscreteModel(parts::AbstractVector{<:Integer},
                                        coarse_model::DiscreteModel{Dc,Dp},
                                        num_uniform_refinements) where {Dc,Dp}
  comm = parts.comm
  if GridapDistributed.i_am_in(comm)
    ptr_pXest_connectivity,
      ptr_pXest,
        ptr_pXest_ghost,
          ptr_pXest_lnodes = setup_ptr_pXest_objects(Val{Dc},
                                                      comm,
                                                      coarse_model,
                                                      num_uniform_refinements)
    dmodel = setup_distributed_discrete_model(Val{Dc},
                                            parts,
                                            coarse_model,
                                            ptr_pXest_connectivity,
                                            ptr_pXest,
                                            ptr_pXest_ghost,
                                            ptr_pXest_lnodes)
    pXest_lnodes_destroy(Val{Dc},ptr_pXest_lnodes)
    pXest_ghost_destroy(Val{Dc},ptr_pXest_ghost)

    non_conforming_glue = _create_conforming_model_non_conforming_glue(dmodel)

    return OctreeDistributedDiscreteModel(Dc,
                                          Dp,
                                          parts,
                                          dmodel,
                                          non_conforming_glue,
                                          coarse_model,
                                          ptr_pXest_connectivity,
                                          ptr_pXest,
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
    pXest_connectivity_destroy(Val{Dc},model.ptr_pXest_connectivity)
  end
  return nothing
end

function octree_distributed_discrete_model_free!(model::OctreeDistributedDiscreteModel{Dc}) where Dc
  if !isa(model.ptr_pXest,Nothing)
    pXest_destroy(Val{Dc},model.ptr_pXest)
  end
  if (model.owns_ptr_pXest_connectivity)
    pXest_connectivity_destroy(Val{Dc},model.ptr_pXest_connectivity)
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

function pXest_copy(::Type{Val{Dc}}, ptr_pXest) where Dc
  if (Dc==2)
    p4est_copy(ptr_pXest, Cint(0))
  else
    p8est_copy(ptr_pXest, Cint(0))
  end
end

function pXest_partition!(::Type{Val{Dc}}, ptr_pXest) where Dc
  if (Dc==2)
    # The 1 here is required to avoid that the children of the 
    # same parent are assigned to different partitions
    p4est_partition(ptr_pXest, 0, C_NULL)
  else
    p8est_partition(ptr_pXest, 0, C_NULL)
  end
end

function pXest_balance!(::Type{Val{Dc}}, ptr_pXest; k_2_1_balance=0) where Dc
  if (Dc==2)
    if (k_2_1_balance==0)
      p4est_balance(ptr_pXest, P4est_wrapper.P4EST_CONNECT_FULL, C_NULL) 
    else 
      p4est_balance(ptr_pXest, P4est_wrapper.P4EST_CONNECT_FACE, C_NULL)
    end
  else
    if (k_2_1_balance==0)
      p8est_balance(ptr_pXest, P4est_wrapper.P8EST_CONNECT_FULL, C_NULL) 
    elseif (k_2_1_balance==1)
      p8est_balance(ptr_pXest, P4est_wrapper.P8EST_CONNECT_EDGE, C_NULL)
    else 
      @assert k_2_1_balance==2
      p8est_balance(ptr_pXest, P4est_wrapper.P8EST_CONNECT_FACE, C_NULL)  
    end
  end
end

function pXest_partition_given!(::Type{Val{Dc}}, ptr_pXest, new_num_cells_per_part) where Dc
  if (Dc==2)
    p4est_partition_given(ptr_pXest, new_num_cells_per_part)
  else
    @assert false
  end
end

function pXest_uniformly_refine!(::Type{Val{Dc}}, ptr_pXest) where Dc
  # Refine callback
  function refine_fn(::Ptr{p4est_t},
                     which_tree::p4est_topidx_t,
                     quadrant::Ptr{p4est_quadrant_t})
    return Cint(1)
  end
  # C-callable refine callback
  refine_fn_c=@cfunction($refine_fn,
                         Cint,
                         (Ptr{p4est_t}, p4est_topidx_t, Ptr{p4est_quadrant_t}))
  if (Dc==2)
    p4est_refine(ptr_pXest, Cint(0), refine_fn_c, C_NULL)
  else
    @assert false
  end
end

function pXest_refine!(::Type{Val{Dc}}, ptr_pXest, refine_fn_c, refine_replace_fn_c; init_fn_c=C_NULL) where Dc
  if (Dc==2)
    p4est_refine_ext(ptr_pXest, Cint(0), Cint(-1), refine_fn_c, init_fn_c, refine_replace_fn_c)
  else
    p8est_refine_ext(ptr_pXest, Cint(0), Cint(-1), refine_fn_c, init_fn_c, refine_replace_fn_c)
  end
end

function pXest_uniformly_coarsen!(::Type{Val{Dc}}, ptr_pXest) where Dc
  # Coarsen callback
  function coarsen_fn(::Ptr{p4est_t},
                      ::p4est_topidx_t,
                      ::Ptr{Ptr{p4est_quadrant_t}})
    return Cint(1)
  end
  # C-callable coasen callback
  coarsen_fn_c=@cfunction($coarsen_fn,
                         Cint,
                         (Ptr{p4est_t}, p4est_topidx_t, Ptr{Ptr{p4est_quadrant_t}}))
  if (Dc==2)
    p4est_coarsen(ptr_pXest, Cint(0), coarsen_fn_c, C_NULL)
  else
    @assert false
  end
end

function get_num_children(::Type{Val{Dc}}) where Dc
  2^Dc
end

function pXest_coarsen!(::Type{Val{Dc}}, ptr_pXest, coarsen_fn_c) where Dc
  if (Dc==2)
    p4est_coarsen(ptr_pXest, Cint(0), coarsen_fn_c, C_NULL)
  else
    p8est_coarsen(ptr_pXest, Cint(0), coarsen_fn_c, C_NULL)
  end
end


function pXest_reset_data!(::Type{Val{Dc}}, ptr_pXest, data_size, init_fn_c, user_pointer) where Dc
  if (Dc==2)
    p4est_reset_data(ptr_pXest, data_size, init_fn_c, user_pointer)
  else
    p8est_reset_data(ptr_pXest, data_size, init_fn_c, user_pointer)
  end
end


# [c/e][child_id][flid]->clid
const rrule_f_to_c_lid_2D=Vector{Vector{Vector{UInt8}}}(
  [
   [[1,1,3,1], [1,2,1,4], [3,1,3,2], [1,4,2,4]],  # c
   [[1,1,3,1], [1,1,1,4], [1,2,3,1], [1,2,1,4]]   # e
  ])

# [c/e][child_id][flid]->clid_dim
const rrule_f_to_c_dim_2D=Vector{Vector{Vector{UInt8}}}(
  [
   [[0,1,1,2], [1,0,2,1], [1,2,0,1], [2,1,1,0]],  # c
   [[1,2,1,2], [1,2,2,1], [2,1,1,2], [2,1,2,1]]   # e
  ])
 
function pXest_update_flags!(::Type{Val{Dc}}, ptr_pXest_old, ptr_pXest_new) where Dc
  pXest_old = ptr_pXest_old[]
  pXest_new = ptr_pXest_new[]
  flags=unsafe_wrap(Array, 
                    Ptr{Cint}(pXest_old.user_pointer), 
                    pXest_old.local_num_quadrants)
  
  num_trees = Cint(pXest_old.connectivity[].num_trees)
  @assert num_trees == Cint(pXest_new.connectivity[].num_trees)

  num_children = get_num_children(Val{Dc})
  global_iquad_new = 0
  global_iquad_old = 0
  for itree = 0:num_trees-1
   tree_old = _pXest_tree_array_index(Val{Dc},pXest_old.trees,itree)[]
   tree_new = _pXest_tree_array_index(Val{Dc},pXest_new.trees,itree)[]
   num_quads_old = Cint(tree_old.quadrants.elem_count)
   local_iquad_old=0
   local_iquad_new=0
   while local_iquad_old < num_quads_old
     q_old = _pXest_quadrant_array_index(Val{Dc},tree_old.quadrants,local_iquad_old)
     q_new = _pXest_quadrant_array_index(Val{Dc},tree_new.quadrants,local_iquad_new)
     if (_pXest_quadrant_compare(Val{Dc},q_old,q_new) == 0)     # q_old was not refined nor coarsened
       flags[global_iquad_old+1] = nothing_flag
       global_iquad_new += 1
       global_iquad_old += 1
       local_iquad_new += 1
       local_iquad_old += 1
     elseif (_pXest_quadrant_is_parent(Val{Dc},q_old,q_new)!=0) # q_old was refined
       flags[global_iquad_old+1] = refine_flag
       global_iquad_new += num_children
       global_iquad_old += 1
       local_iquad_new += num_children
       local_iquad_old += 1
     elseif (_pXest_quadrant_is_parent(Val{Dc},q_new,q_old)!=0) # q_old and its siblings were coarsened 
       for i=0:num_children-1
         flags[global_iquad_old+i+1] = coarsen_flag
       end
       global_iquad_old += num_children
       global_iquad_new += 1
       local_iquad_old += num_children
       local_iquad_new += 1
     else
       @assert false
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
  f1,f2,f3, cgids_snd, cgids_rcv = map(fmodel.models,
                                       partition(fgids)) do fmodel, fpartition
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
  cache=fetch_vector_ghost_values_cache(dfcell_to_child_id,partition(fgids))
  fetch_vector_ghost_values!(dfcell_to_child_id,cache) |> wait
  
  # Note: Reversing snd and rcv 
  parts_rcv, parts_snd = assembly_neighbors(partition(fgids))
  PArrays.exchange_fetch!(cgids_rcv,cgids_snd,ExchangeGraph(parts_snd,parts_rcv))

  map(f1,cgids_rcv) do fine_to_coarse_faces_map, cgids_rcv
    if (GridapDistributed.i_am_in(cparts))
      cgids=get_cell_gids(cmodel)
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
  map(f1,f2,f3) do fine_to_coarse_faces_map, fine_to_coarse_faces_dim, fcell_to_child_id
    if (!(GridapDistributed.i_am_in(cparts)))
      nothing
    else
      polytope=(Dc==2 ? QUAD : HEX)
      cmodel_local = PArrays.getany(cmodel.models)
      num_cells_coarse = num_cells(cmodel_local)
      reffe  = LagrangianRefFE(Float64,polytope,1)
      rrules = Fill(Gridap.Adaptivity.RefinementRule(reffe,2),num_cells_coarse)
      AdaptivityGlue(fine_to_coarse_faces_map,fcell_to_child_id,rrules)
    end
  end
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

function _move_fwd_and_check_if_all_children_coarsened(flags,num_o_c_cells,cell,num_children)
  e=cell+num_children-1
  while (cell <= num_o_c_cells) && (cell <= e)
    if (flags[cell]!=coarsen_flag)
      break
    end
    cell=cell+1
  end
  return cell,cell==e+1
end

function _compute_fine_to_coarse_model_glue(
         cparts,
         cmodel::Union{Nothing,GridapDistributed.DistributedDiscreteModel{Dc}},
         fmodel::GridapDistributed.DistributedDiscreteModel{Dc},
         refinement_and_coarsening_flags::MPIArray{<:AbstractVector}) where Dc

  function setup_communication_buffers_fine_partition(cparts,
                                              fmodel,
                                              cmodel::Union{Nothing,GridapDistributed.DistributedDiscreteModel},
                                              fine_to_coarse_faces_map::Union{Nothing,MPIArray},
                                              fcell_to_child_id::Union{Nothing,MPIArray})
    fgids=get_cell_gids(fmodel)
    map(fmodel.models,fine_to_coarse_faces_map,fcell_to_child_id) do fmodel, fine_to_coarse_faces_map, fcell_to_child_id
      if (!(GridapDistributed.i_am_in(cparts)))
        # cmodel might be distributed among less processes than fmodel
        Int[], Int[], Int[], Int[]
      else
        cgids        = get_cell_gids(cmodel)
        cpartition   = PArrays.getany(partition(cgids))
        _setup_communication_buffers_fine_partition(fine_to_coarse_faces_map,
                                                    fcell_to_child_id,
                                                    partition(fgids),
                                                    cpartition)
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
    cgids_data  = local_to_global(cpartition)[fine_to_coarse_faces_map[end][lids_snd.data]]
    cgids_snd   = PArrays.JaggedArray(cgids_data,lids_snd.ptrs)
    cgids_rcv   = PArrays.JaggedArray(Vector{Int}(undef,length(lids_rcv.data)),lids_rcv.ptrs)
    cell_child_id_data = fcell_to_child_id[lids_snd.data]
    cell_child_id_snd = PArrays.JaggedArray(cell_child_id_data,lids_snd.ptrs)
    fcell_child_id_rcv = PArrays.JaggedArray(Vector{Int}(undef,length(lids_rcv.data)),lids_rcv.ptrs)
    cgids_snd, cgids_rcv, cell_child_id_snd, fcell_child_id_rcv
  end

  function _setup_communication_buffers_fine_partition(
    fine_to_coarse_faces_map::AbstractVector{<:Gridap.Arrays.Table},
    fcell_to_child_id::Gridap.Arrays.Table,
    fpartition::MPIArray,
    cpartition::AbstractLocalIndices)

    # Note: Reversing snd and rcv
    lids_rcv,lids_snd = map(PArrays.getany,assembly_local_indices(fpartition))

    cgids_data   = Vector{Int}(undef,length(lids_snd.data))
    cgids_data  .= -1 

    fcell_child_id_data = Vector{Int}(undef,length(lids_snd.data))
    fcell_child_id_data .= -1
    
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
    cgids_snd   = PArrays.JaggedArray(cgids_data,lids_snd.ptrs)
    cgids_rcv   = PArrays.JaggedArray(Vector{Int}(undef,length(lids_rcv.data)),lids_rcv.ptrs)
    cell_child_id_snd = PArrays.JaggedArray(fcell_child_id_data,lids_snd.ptrs)
    cell_child_id_rcv = PArrays.JaggedArray(Vector{Int}(undef,length(lids_rcv.data)),lids_rcv.ptrs)
    cgids_snd, cgids_rcv, cell_child_id_snd, cell_child_id_rcv
  end

  function setup_communication_buffers_coarse_partition(cparts,
                                              cmodel::Union{Nothing,GridapDistributed.DistributedDiscreteModel},
                                              fmodel,
                                              fine_to_coarse_faces_map::Union{Nothing,MPIArray},
                                              fcell_to_child_id::Union{Nothing,MPIArray})
    fgids=get_cell_gids(fmodel)
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

    fgids_data   = Vector{Int}(undef,length(clids_snd.data))
    fgids_data  .= -1
    
    child_ids_data = Vector{Int}(undef,length(clids_snd.data))
    child_ids_data .= -1

    f2c_map_ptrs = fine_to_coarse_faces_map[end].ptrs 
    f2c_map_data = fine_to_coarse_faces_map[end].data
    fchild_id_data = fcell_to_child_id.data
    fl2g         = local_to_global(PArrays.getany(fpartition))

    for i=1:length(flids_snd.ptrs)-1
      for j=flids_snd.ptrs[i]:flids_snd.ptrs[i+1]-1
        fcell = flids_snd.data[j]
        # fcell coarsened ...
        if (f2c_map_ptrs[fcell+1]-f2c_map_ptrs[fcell]>1)
          for k=f2c_map_ptrs[fcell]:f2c_map_ptrs[fcell+1]-1
              ccell = f2c_map_data[k]
              child_id = fchild_id_data[k]
              # find ccell in clids_snd[i]:clids_snd[i+1]-1
              for l=clids_snd.ptrs[i]:clids_snd.ptrs[i+1]-1
                if (clids_snd.data[l]==ccell)
                  fgids_data[l]  = fl2g[fcell]
                  child_ids_data[l] = child_id 
                end
              end
          end
        end   
      end 
    end
    fgids_snd     = PArrays.JaggedArray(fgids_data,clids_snd.ptrs)
    fgids_rcv     = PArrays.JaggedArray(Vector{Int}(undef,length(clids_rcv.data)),clids_rcv.ptrs)
    child_id_snd = PArrays.JaggedArray(child_ids_data,clids_snd.ptrs)
    child_id_rcv = PArrays.JaggedArray(Vector{Int}(undef,length(clids_rcv.data)),clids_rcv.ptrs)
    fgids_snd, fgids_rcv, child_id_snd, child_id_rcv
  end

  function _check_if_coarsen(Dc,
                            cpartition,
                            flags)
    num_children = get_num_children(Val{Dc})
    num_o_c_cells = own_length(cpartition)
    cell = 1
    while cell <= num_o_c_cells
      if (flags[cell]==coarsen_flag)
        cell,coarsen=_move_fwd_and_check_if_all_children_coarsened(flags,num_o_c_cells,cell,num_children)
        if coarsen
          return true
        end 
      else 
        cell=cell+1
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

          for i=1:length(flids_snd.ptrs)-1
            for j=flids_snd.ptrs[i]:flids_snd.ptrs[i+1]-1
              fcell = flids_snd.data[j]
              # fcell coarsened ...
              if (f2c_map_ptrs[fcell+1]-f2c_map_ptrs[fcell]>1)
                for k=f2c_map_ptrs[fcell]:f2c_map_ptrs[fcell+1]-1
                   ccell = f2c_map_data[k]
                   # find ccell in clids_snd[i]:clids_snd[i+1]-1
                   for l=clids_snd.ptrs[i]:clids_snd.ptrs[i+1]-1
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
          snd_buffer   = PArrays.JaggedArray(snd_buffer_data,flids_snd.ptrs)
          rcv_buffer   = PArrays.JaggedArray(Vector{Int}(undef,length(flids_rcv.data)),flids_rcv.ptrs)
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
          cgids=get_cell_gids(cmodel)
          cpartition = PArrays.getany(partition(cgids))
          # Note: Reversing snd and rcv
          lids_rcv,_ = map(PArrays.getany,assembly_local_indices(partition(fgids)))
          glo_to_loc = global_to_local(cpartition)
          for i in 1:length(lids_rcv.data)
            lid = lids_rcv.data[i]
            gid = cgids_rcv.data[i]
            fine_to_coarse_faces_map[end][lid] = glo_to_loc[gid]
            fcell_to_child_id[lid] = child_id_rcv.data[i]
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
          data=fine_to_coarse_faces_map[end].data
          ptrs=fine_to_coarse_faces_map[end].ptrs
          cgids=get_cell_gids(cmodel)
          cpartition = PArrays.getany(partition(cgids))
          # Note: Reversing snd and rcv
          lids_rcv,_ = map(PArrays.getany,assembly_local_indices(partition(fgids)))
          glo_to_loc = global_to_local(cpartition)
          for i in 1:length(lids_rcv.data)
            lid = lids_rcv.data[i]
            gid = cgids_rcv.data[i]
            if gid!=-1
               @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))] gid=$(gid) $(ptrs[lid+1]-ptrs[lid])"
               @assert (ptrs[lid+1]-ptrs[lid])==1
               data[ptrs[lid]] = glo_to_loc[gid]
               fcell_to_child_id.data[ptrs[lid]] = child_id_rcv.data[i]
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
          f2c_data=fine_to_coarse_faces_map[end].data
          f2c_ptrs=fine_to_coarse_faces_map[end].ptrs
          fchild_id_data=fcell_to_child_id.data
          cgids=get_cell_gids(cmodel)
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
  coarsen_array=map(partition(cgids),refinement_and_coarsening_flags) do cpartition,flags
    _check_if_coarsen(Dc,cpartition,flags) 
  end 
  or_func(a,b)=a || b
  coarsen=PArrays.getany(reduction(or_func,coarsen_array;destination=:all,init=false))

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
          _process_owned_cells_fine_to_coarse_model_glue(cmodel_local,fmodel,cpartition,fpartition,flags,coarsen)
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
  cache=fetch_vector_ghost_values_cache(refinement_and_coarsening_flags,partition(cgids))
  fetch_vector_ghost_values!(refinement_and_coarsening_flags,cache) |> wait
  
  # Note: Reversing snd and rcv 
  fparts_rcv, fparts_snd = assembly_neighbors(partition(fgids))
  graph=ExchangeGraph(fparts_snd,fparts_rcv)
  PArrays.exchange_fetch!(cgids_rcv,cgids_snd,graph)
  PArrays.exchange_fetch!(fchild_id_rcv,fchild_id_snd,graph)

  if coarsen 
    cparts_rcv, cparts_snd = assembly_neighbors(partition(cgids))
    graph=ExchangeGraph(cparts_snd,cparts_rcv)
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
  map(f1,f2,f3,refinement_and_coarsening_flags) do fine_to_coarse_faces_map, 
                                                   fine_to_coarse_faces_dim, 
                                                   fcell_to_child_id,
                                                   flags
    if (!(GridapDistributed.i_am_in(cparts)))
      nothing
    else
      ## The following lines are a replacement for WhiteRefinementRule()
      ## to have the types of rrule_nothing_flag and rrule_refinement_flag
      ## to be 100% equivalent for all type parameters
      polytope=(Dc==2 ? QUAD : HEX)
      partition = Gridap.ReferenceFEs.tfill(1,Val{Dc}())
      ref_grid = UnstructuredGrid(compute_reference_grid(polytope,partition))
      rrule_nothing_flag = 
         Gridap.Adaptivity.RefinementRule(Gridap.Adaptivity.WithoutRefinement(),polytope,ref_grid)
      
      reffe  = LagrangianRefFE(Float64,polytope,1)
      rrule_refinement_flag = 
          Gridap.Adaptivity.RefinementRule(reffe,2)
      f(x)=x==nothing_flag ? 1 : 2
      coarse_cell_to_rrule=map(f,flags)
      rrules=Gridap.Arrays.CompressedArray([rrule_nothing_flag,rrule_refinement_flag],coarse_cell_to_rrule)
      @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))] fine_to_coarse_faces_map[end]: $(fine_to_coarse_faces_map[end])"
      @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))] fcell_to_child_id: $(fcell_to_child_id)"
      GT=isa(fine_to_coarse_faces_map,Vector{<:AbstractVector{<:Integer}}) ? Gridap.Adaptivity.RefinementGlue() : Gridap.Adaptivity.MixedGlue()
      AdaptivityGlue(GT,fine_to_coarse_faces_map,fcell_to_child_id,rrules)
    end
  end
end

function _process_owned_cells_fine_to_coarse_model_glue(cmodel::DiscreteModel{Dc},
                                                        fmodel::DiscreteModel{Dc},
                                                        cpartition,
                                                        fpartition,
                                                        flags,
                                                        coarsen) where Dc

  function _setup_fine_to_coarse_faces_map_table(Dc,flags,num_o_c_cells,num_f_cells)
    num_children = get_num_children(Val{Dc})
    fine_to_coarse_faces_map_ptrs = Vector{Int}(undef,num_f_cells+1)
    fine_to_coarse_faces_map_ptrs[1]=1
    cell=1
    c = 1
    while cell <= num_o_c_cells
      if (flags[cell]==nothing_flag) 
        fine_to_coarse_faces_map_ptrs[c+1]=fine_to_coarse_faces_map_ptrs[c]+1
        cell=cell+1
        c=c+1
      elseif (flags[cell]==refine_flag)
        for child=1:num_children
          fine_to_coarse_faces_map_ptrs[c+1]=fine_to_coarse_faces_map_ptrs[c]+1
          c=c+1
        end
        cell=cell+1
      else
        @assert flags[cell]==coarsen_flag
        cell_fwd,coarsen=_move_fwd_and_check_if_all_children_coarsened(flags,num_o_c_cells,cell,num_children)
        if coarsen
          @assert cell_fwd-cell==num_children
          fine_to_coarse_faces_map_ptrs[c+1]=fine_to_coarse_faces_map_ptrs[c]+num_children
          cell=cell+num_children
          c=c+1
        else 
          for j=c:c+(cell_fwd-cell+1) 
            fine_to_coarse_faces_map_ptrs[j+1]=fine_to_coarse_faces_map_ptrs[j]+1
            c=c+1
          end
          cell=cell_fwd
          c=c+cell_fwd-cell
        end
      end
    end
    for j=c:num_f_cells
      fine_to_coarse_faces_map_ptrs[j+1]=fine_to_coarse_faces_map_ptrs[j]
    end
    fine_to_coarse_faces_map_data = Vector{Int}(undef,fine_to_coarse_faces_map_ptrs[end]-1)
    fcell_to_child_id_data = Vector{Int}(undef,fine_to_coarse_faces_map_ptrs[end]-1)
    cell=1
    c = 1
    while cell <= num_o_c_cells
      if (flags[cell]==refine_flag)
        for child = 1:num_children
          fine_to_coarse_faces_map_data[fine_to_coarse_faces_map_ptrs[c]+child-1] = cell
          fcell_to_child_id_data[fine_to_coarse_faces_map_ptrs[c]+child-1] = child
        end 
        c = c + num_children
        cell=cell+1
      elseif (flags[cell]==nothing_flag)
        fine_to_coarse_faces_map_data[fine_to_coarse_faces_map_ptrs[c]]=cell
        fcell_to_child_id_data[fine_to_coarse_faces_map_ptrs[c]]=1
        c=c+1
        cell=cell+1
      else
        @assert flags[cell]==coarsen_flag
        cell_fwd,coarsen=_move_fwd_and_check_if_all_children_coarsened(flags,num_o_c_cells,cell,num_children)
        if coarsen
          for child = 1:num_children
            fcell_to_child_id_data[fine_to_coarse_faces_map_ptrs[c]+child-1] = child
            fine_to_coarse_faces_map_data[fine_to_coarse_faces_map_ptrs[c]+child-1] = cell
            cell=cell+1
          end 
          c=c+1
        else 
          for j=cell:cell_fwd-1 
            fine_to_coarse_faces_map_data[fine_to_coarse_faces_map_ptrs[c]]=j
            fcell_to_child_id_data[fine_to_coarse_faces_map_ptrs[c]]=1
            c=c+1
            cell=cell+1
          end
        end
      end
    end
    Gridap.Arrays.Table(fine_to_coarse_faces_map_data, fine_to_coarse_faces_map_ptrs),
             Gridap.Arrays.Table(fcell_to_child_id_data, fine_to_coarse_faces_map_ptrs)
  end

  function _setup_fine_to_coarse_faces_map_vector!(fine_to_coarse_faces_map,fcell_to_child_id,Dc,flags,num_o_c_cells)
    # Go over all cells of coarse grid portion
    num_children = get_num_children(Val{Dc})
    c = 1
    for cell = 1:num_o_c_cells
      if flags[cell]==refine_flag
        for child = 1:num_children
          fine_to_coarse_faces_map[c+child-1] = cell
          fcell_to_child_id[c+child-1] = child
        end
        c = c + num_children
      elseif (flags[cell]==nothing_flag)
        fine_to_coarse_faces_map[c] = cell
        fcell_to_child_id[c] = 1
        c=c+1
      else
        @assert flags[cell]!=coarsen_flag 
        error("Unknown AMR flag")
      end
    end
  end
  
  num_f_cells   = num_cells(fmodel)       # Number of fine cells (owned+ghost)
  num_o_c_cells = own_length(cpartition)  # Number of coarse cells (owned)
  ftopology = Gridap.Geometry.get_grid_topology(fmodel)
  if coarsen
    fine_to_coarse_faces_map = Vector{Gridap.Arrays.Table{Int,Vector{Int},Vector{Int}}}(undef,Dc+1)
    a,b=_setup_fine_to_coarse_faces_map_table(Dc,flags,num_o_c_cells,num_f_cells)
    fine_to_coarse_faces_map[Dc+1]=a
    fcell_to_child_id=b
    # In the future we should also have here the code to also setup
    # fine_to_coarse_faces_map for lower dimensional objects
  else
    fine_to_coarse_faces_map = Vector{Vector{Int}}(undef,Dc+1)
    fine_to_coarse_faces_map[Dc+1] = Vector{Int}(undef,num_f_cells)
    fcell_to_child_id = Vector{Int}(undef,num_f_cells)
    _setup_fine_to_coarse_faces_map_vector!(fine_to_coarse_faces_map[Dc+1],fcell_to_child_id,Dc,flags,num_o_c_cells)
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
     ptr_new_pXest = pXest_copy(Val{Dc}, model.ptr_pXest)
     pXest_uniformly_refine!(Val{Dc}, ptr_new_pXest)
   else
     ptr_new_pXest = nothing
   end

   new_comm = isa(parts,Nothing) ? old_comm : parts.comm
   if GridapDistributed.i_am_in(new_comm)
      if !isa(parts,Nothing)
        aux = ptr_new_pXest
        ptr_new_pXest = _p4est_to_new_comm(ptr_new_pXest,
                                           model.ptr_pXest_connectivity,
                                           model.parts.comm,
                                           parts.comm)
        if GridapDistributed.i_am_in(old_comm)
          pXest_destroy(Val{Dc},aux)
        end
      end

      # Extract ghost and lnodes
      ptr_pXest_ghost  = setup_pXest_ghost(Val{Dc}, ptr_new_pXest)
      ptr_pXest_lnodes = setup_pXest_lnodes(Val{Dc}, ptr_new_pXest, ptr_pXest_ghost)

      # Build fine-grid mesh
      new_parts = isa(parts,Nothing) ? model.parts : parts
      fmodel = setup_distributed_discrete_model(Val{Dc},
                                              new_parts,
                                              model.coarse_model,
                                              model.ptr_pXest_connectivity,
                                              ptr_new_pXest,
                                              ptr_pXest_ghost,
                                              ptr_pXest_lnodes)

      pXest_lnodes_destroy(Val{Dc},ptr_pXest_lnodes)
      pXest_ghost_destroy(Val{Dc},ptr_pXest_ghost)

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
                                     false,
                                     model)


      return ref_model, dglue
   else
    new_parts = isa(parts,Nothing) ? model.parts : parts
    return VoidOctreeDistributedDiscreteModel(model,new_parts), nothing
   end
end

function Gridap.Adaptivity.adapt(model::OctreeDistributedDiscreteModel{Dc,Dp}, 
                                  refinement_and_coarsening_flags::MPIArray{<:Vector};
                                  parts=nothing) where {Dc,Dp}

    Gridap.Helpers.@notimplementedif parts!=nothing

    # Variables which are updated accross calls to init_fn_callback_2d
    current_quadrant_index_within_tree = Cint(0)
    current_quadrant_index_among_trees = Cint(0)

    # This C callback function is called once per quadtree quadrant. Here we are assuming
    # that p4est->user_pointer has been set prior to the first call to this call
    # back function to an array of ints with as many entries as forest quadrants. This call back function
    # initializes the quadrant->p.user_data void * pointer of all quadrants such that it
    # points to the corresponding entry in the global array mentioned in the previous sentence.
    if (Dc==2)
      function init_fn_callback_2d(forest_ptr::Ptr{p4est_t},
        which_tree::p4est_topidx_t,
        quadrant_ptr::Ptr{p4est_quadrant_t})
        # Extract a reference to the tree which_tree
        forest = forest_ptr[]
        tree = p4est_tree_array_index(forest.trees, which_tree)[]
        quadrant = quadrant_ptr[]
        q = P4est_wrapper.p4est_quadrant_array_index(tree.quadrants, current_quadrant_index_within_tree)
        @assert p4est_quadrant_compare(q, quadrant_ptr) == 0
        user_data = unsafe_wrap(Array, 
                                Ptr{Cint}(forest.user_pointer), 
                                current_quadrant_index_among_trees+1)[current_quadrant_index_among_trees+1]
        unsafe_store!(Ptr{Cint}(quadrant.p.user_data), user_data, 1)
        current_quadrant_index_within_tree = (current_quadrant_index_within_tree + 1) % (tree.quadrants.elem_count)
        current_quadrant_index_among_trees = current_quadrant_index_among_trees+1
        return nothing
      end
      init_fn_callback_2d_c = @cfunction($init_fn_callback_2d, 
                                        Cvoid, (Ptr{p4est_t}, p4est_topidx_t, Ptr{p4est_quadrant_t}))
      init_fn_callback_c = init_fn_callback_2d_c


      function coarsen_callback_2d(forest_ptr::Ptr{p4est_t},
                                   which_tree::p4est_topidx_t,
                                   quadrant_ptr::Ptr{Ptr{p4est_quadrant_t}})

        num_children=get_num_children(Val{2})
        quadrants=unsafe_wrap(Array, quadrant_ptr, num_children)
        coarsen=Cint(1)
        for quadrant_index=1:num_children
          quadrant = quadrants[quadrant_index][]
          # I have noticed that new quadrants created as by-product
          # of the refininement process have quadrant.p.user_data == C_NULL
          # Not sure why ... The following if-end takes care of this.
          if (quadrant.p.user_data) == C_NULL
            return Cint(0)
          end
          is_coarsen_flag=(unsafe_wrap(Array,Ptr{Cint}(quadrant.p.user_data),1)[])==coarsen_flag
          if (!is_coarsen_flag) 
            return Cint(0)
          end
        end  
        return coarsen
      end 
      coarsen_fn_callback_2d_c = @cfunction($coarsen_callback_2d, 
                                            Cint, (Ptr{p4est_t}, p4est_topidx_t, Ptr{Ptr{p4est_quadrant_t}}))
      coarsen_fn_callback_c = coarsen_fn_callback_2d_c


      function refine_replace_callback_2d(::Ptr{p4est_t},
                                        which_tree::p4est_topidx_t,
                                        num_outgoing::Cint,
                                        outgoing_ptr::Ptr{Ptr{p4est_quadrant_t}},
                                        num_incoming::Cint,
                                        incoming_ptr::Ptr{Ptr{p4est_quadrant_t}})
        num_children=get_num_children(Val{2}) 
        @assert num_outgoing==1 
        @assert num_incoming==num_children
        outgoing=unsafe_wrap(Array, outgoing_ptr, 1)
        quadrant = outgoing[1][]
        incoming=unsafe_wrap(Array, incoming_ptr, num_children)
        for quadrant_index=1:num_children
          quadrant = incoming[quadrant_index][]
          if (quadrant.p.user_data) != C_NULL
             unsafe_store!(Ptr{Cint}(quadrant.p.user_data), nothing_flag, 1)
          end
        end
     end

    refine_replace_callback_2d_c = 
       @cfunction($refine_replace_callback_2d, Cvoid, (Ptr{p4est_t}, 
                                                       p4est_topidx_t, 
                                                       Cint, 
                                                       Ptr{Ptr{p4est_quadrant_t}}, 
                                                       Cint, 
                                                       Ptr{Ptr{p4est_quadrant_t}}))
    
    refine_replace_callback_c = refine_replace_callback_2d_c     

    else
      @assert Dc==3
      function init_fn_callback_3d(forest_ptr::Ptr{p8est_t},
        which_tree::p4est_topidx_t,
        quadrant_ptr::Ptr{p8est_quadrant_t})
        # Extract a reference to the tree which_tree
        forest = forest_ptr[]
        tree = p8est_tree_array_index(forest.trees, which_tree)[]
        quadrant = quadrant_ptr[]
        q = P4est_wrapper.p8est_quadrant_array_index(tree.quadrants, current_quadrant_index_within_tree)
        @assert p8est_quadrant_compare(q, quadrant_ptr) == 0
        user_data = unsafe_wrap(Array, 
                                Ptr{Cint}(forest.user_pointer), 
                                current_quadrant_index_among_trees+1)[current_quadrant_index_among_trees+1]
        unsafe_store!(Ptr{Cint}(quadrant.p.user_data), user_data, 1)
        current_quadrant_index_within_tree = (current_quadrant_index_within_tree + 1) % (tree.quadrants.elem_count)
        current_quadrant_index_among_trees = current_quadrant_index_among_trees+1
        return nothing
      end
      init_fn_callback_3d_c = @cfunction($init_fn_callback_3d, 
                                        Cvoid, (Ptr{p8est_t}, p4est_topidx_t, Ptr{p8est_quadrant_t}))
      init_fn_callback_c = init_fn_callback_3d_c


       function coarsen_callback_3d(forest_ptr::Ptr{p8est_t},
                                   which_tree::p4est_topidx_t,
                                   quadrant_ptr::Ptr{Ptr{p8est_quadrant_t}})

        num_children=get_num_children(Val{3})
        quadrants=unsafe_wrap(Array, quadrant_ptr, num_children)
        coarsen=Cint(1)
        for quadrant_index=1:num_children
          quadrant = quadrants[quadrant_index][]
          # I have noticed that new quadrants created as by-product
          # of the refininement process have quadrant.p.user_data == C_NULL
          # Not sure why ... The following if-end takes care of this.
          if (quadrant.p.user_data) == C_NULL
            return Cint(0)
          end
          is_coarsen_flag=(unsafe_wrap(Array,Ptr{Cint}(quadrant.p.user_data),1)[])==coarsen_flag
          if (!is_coarsen_flag) 
            return Cint(0)
          end
        end  
        return coarsen
      end 
      coarsen_fn_callback_3d_c = @cfunction($coarsen_callback_3d, 
                                            Cint, (Ptr{p8est_t}, p4est_topidx_t, Ptr{Ptr{p8est_quadrant_t}}))
      coarsen_fn_callback_c = coarsen_fn_callback_3d_c

      function refine_replace_callback_3d(::Ptr{p8est_t},
                                        which_tree::p4est_topidx_t,
                                        num_outgoing::Cint,
                                        outgoing_ptr::Ptr{Ptr{p8est_quadrant_t}},
                                        num_incoming::Cint,
                                        incoming_ptr::Ptr{Ptr{p8est_quadrant_t}})
        num_children=get_num_children(Val{3}) 
        @assert num_outgoing==1 
        @assert num_incoming==num_children
        outgoing=unsafe_wrap(Array, outgoing_ptr, 1)
        quadrant = outgoing[1][]
        incoming=unsafe_wrap(Array, incoming_ptr, num_children)
        for quadrant_index=1:num_children
          quadrant = incoming[quadrant_index][]
          if (quadrant.p.user_data) != C_NULL
             unsafe_store!(Ptr{Cint}(quadrant.p.user_data), nothing_flag, 1)
          end
        end
     end

     refine_replace_callback_3d_c = 
       @cfunction($refine_replace_callback_3d, Cvoid, (Ptr{p8est_t}, 
                                                       p4est_topidx_t, 
                                                       Cint, 
                                                       Ptr{Ptr{p8est_quadrant_t}}, 
                                                       Cint, 
                                                       Ptr{Ptr{p8est_quadrant_t}}))

     refine_replace_callback_c = refine_replace_callback_3d_c
    end                                     

    map(model.dmodel.models,refinement_and_coarsening_flags) do lmodel, flags
      # The length of the local flags array has to match the number of 
      # cells in the model. This includes both owned and ghost cells. 
      # Only the flags for owned cells are actually taken into account. 
      @assert num_cells(lmodel)==length(flags)
      pXest_reset_data!(Val{Dc}, model.ptr_pXest, Cint(sizeof(Cint)), init_fn_callback_c, pointer(flags))
    end
    
    function refine_callback_2d(::Ptr{p4est_t},
      which_tree::p4est_topidx_t,
      quadrant_ptr::Ptr{p4est_quadrant_t})
      quadrant = quadrant_ptr[]
      return Cint(unsafe_wrap(Array, Ptr{Cint}(quadrant.p.user_data), 1)[] == refine_flag)
    end
    refine_callback_2d_c = @cfunction($refine_callback_2d, Cint, (Ptr{p4est_t}, p4est_topidx_t, Ptr{p4est_quadrant_t}))



    
    # Copy input p4est, refine and balance
    ptr_new_pXest = pXest_copy(Val{Dc}, model.ptr_pXest)
    pXest_refine!(Val{Dc}, ptr_new_pXest,
                  refine_callback_2d_c,
                  refine_replace_callback_c)
    pXest_coarsen!(Val{Dc}, ptr_new_pXest, coarsen_fn_callback_c)
    pXest_balance!(Val{Dc}, ptr_new_pXest)
    pXest_update_flags!(Val{Dc},model.ptr_pXest,ptr_new_pXest)

    # Extract ghost and lnodes
    ptr_pXest_ghost  = setup_pXest_ghost(Val{Dc}, ptr_new_pXest)
    ptr_pXest_lnodes = setup_pXest_lnodes_nonconforming(Val{Dc}, ptr_new_pXest, ptr_pXest_ghost)

    # Build fine-grid mesh
    fmodel,non_conforming_glue=setup_non_conforming_distributed_discrete_model(Val{Dc},
                                                                               model.parts,
                                                                               model.coarse_model,
                                                                               model.ptr_pXest_connectivity,
                                                                               ptr_new_pXest,
                                                                               ptr_pXest_ghost,
                                                                               ptr_pXest_lnodes)
     
     pXest_ghost_destroy(Val{Dc},ptr_pXest_ghost)
     pXest_lnodes_destroy(Val{Dc},ptr_pXest_lnodes)

     adaptivity_glue = _compute_fine_to_coarse_model_glue(model.parts,
                                                          model.dmodel,
                                                          fmodel,
                                                          refinement_and_coarsening_flags)
     adaptive_models=map(local_views(model),
                         local_views(fmodel),
                         adaptivity_glue) do model, fmodel, glue 
        Gridap.Adaptivity.AdaptedDiscreteModel(fmodel,model,glue)
     end
     fmodel=GridapDistributed.GenericDistributedDiscreteModel(adaptive_models,get_cell_gids(fmodel))

     ref_model = OctreeDistributedDiscreteModel(Dc,Dp,
                                    model.parts,
                                    fmodel,
                                    non_conforming_glue,
                                    model.coarse_model,
                                    model.ptr_pXest_connectivity,
                                    ptr_new_pXest,
                                    false,
                                    model)
     return ref_model, adaptivity_glue
  # else
  #  new_parts = isa(parts,Nothing) ? model.parts : parts
  #  return VoidOctreeDistributedDiscreteModel(model,new_parts), nothing
  # end
end

function Gridap.Adaptivity.coarsen(model::OctreeDistributedDiscreteModel{Dc,Dp}) where {Dc,Dp}
  comm = model.parts.comm
  if (GridapDistributed.i_am_in(comm))
    # Copy and refine input p4est
    ptr_new_pXest = pXest_copy(Val{Dc}, model.ptr_pXest)
    pXest_uniformly_coarsen!(Val{Dc}, ptr_new_pXest)
  else
    ptr_new_pXest=nothing
  end

  if (GridapDistributed.i_am_in(comm))
     # Extract ghost and lnodes
     ptr_pXest_ghost  = setup_pXest_ghost(Val{Dc}, ptr_new_pXest)
     ptr_pXest_lnodes = setup_pXest_lnodes(Val{Dc}, ptr_new_pXest, ptr_pXest_ghost)

     # Build coarse-grid mesh
     cmodel = setup_distributed_discrete_model(Val{Dc},
                                             model.parts,
                                             model.coarse_model,
                                             model.ptr_pXest_connectivity,
                                             ptr_new_pXest,
                                             ptr_pXest_ghost,
                                             ptr_pXest_lnodes)

     pXest_lnodes_destroy(Val{Dc},ptr_pXest_lnodes)
     pXest_ghost_destroy(Val{Dc},ptr_pXest_ghost)

     dglue = _compute_fine_to_coarse_model_glue(model.parts,
                                                cmodel,
                                                model.dmodel)

     nc_glue=_create_conforming_model_non_conforming_glue(cmodel)

     c_octree_model = OctreeDistributedDiscreteModel(Dc,Dp,
                                    model.parts,
                                    cmodel,
                                    nc_glue,
                                    model.coarse_model,
                                    model.ptr_pXest_connectivity,
                                    ptr_new_pXest,
                                    false,
                                    model)
     return c_octree_model, dglue
  else
     return VoidOctreeDistributedDiscreteModel(model,model.parts), nothing
  end
end


# We have a p4est distributed among P processors. This function
# instantiates the same among Q processors.
function _p4est_to_new_comm(ptr_pXest, ptr_pXest_conn, old_comm, new_comm)
  A=is_included(old_comm,new_comm) # old \subset new (smaller to larger nparts)
  B=is_included(new_comm,old_comm) # old \supset new (larger to smaller nparts)
  @assert xor(A,B)
  if (A)
    _p4est_to_new_comm_old_subset_new(ptr_pXest, ptr_pXest_conn, old_comm, new_comm)
  else
    _p4est_to_new_comm_old_supset_new(ptr_pXest, ptr_pXest_conn, old_comm, new_comm)
  end
end

function _p4est_to_new_comm_old_subset_new(ptr_pXest, ptr_pXest_conn, old_comm, new_comm)
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
      quadrants = P4est_wrapper.p4est_deflate_quadrants(ptr_pXest,C_NULL)
      p4est_comm_count_pertree(ptr_pXest,pertree)
      MPI.Bcast!(pertree,0,new_comm)
    else
      MPI.Bcast!(global_first_quadrant,0,new_comm)
      quadrants = sc_array_new_count(sizeof(p4est_quadrant_t), 0)
      MPI.Bcast!(pertree,0,new_comm)
    end
    return P4est_wrapper.p4est_inflate(new_comm,
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

function _p4est_to_new_comm_old_supset_new(ptr_pXest, ptr_pXest_conn, old_comm, new_comm)
  @assert GridapDistributed.i_am_in(old_comm)
  pXest = ptr_pXest[]
  pXest_conn = ptr_pXest_conn[]

  pertree = Vector{P4est_wrapper.p4est_gloidx_t}(undef,pXest_conn.num_trees+1)
  p4est_comm_count_pertree(ptr_pXest,pertree)

  if (GridapDistributed.i_am_in(new_comm))
    new_comm_num_parts = GridapDistributed.num_parts(new_comm)
    global_first_quadrant = Vector{P4est_wrapper.p4est_gloidx_t}(undef,new_comm_num_parts+1)
    pXest=ptr_pXest[]
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
    quadrants = P4est_wrapper.p4est_deflate_quadrants(ptr_pXest,C_NULL)

    return P4est_wrapper.p4est_inflate(new_comm,
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


function _pXest_tree_array_index(::Type{Val{Dc}},trees,itree) where Dc
  if (Dc==2)
    return p4est_tree_array_index(trees,itree)
  elseif (Dc==3)
    return p8est_tree_array_index(trees,itree)
  end
end

function _pXest_comm_find_owner(::Type{Val{Dc}},ptr_pXest,itree,quad,guess) where Dc
  if (Dc==2)
    return p4est_comm_find_owner(ptr_pXest,itree,quad,guess)
  elseif (Dc==3)
    return p8est_comm_find_owner(ptr_pXest,itree,quad,guess)
  end
end

function _pXest_quadrant_array_index(::Type{Val{Dc}}, quadrants, iquad) where Dc
  if (Dc==2)
    return p4est_quadrant_array_index(quadrants, iquad)
  elseif (Dc==3)
    return p8est_quadrant_array_index(quadrants, iquad)
  end
end


function _pXest_quadrant_is_equal(::Type{Val{Dc}},  q1, q2) where Dc
  if (Dc==2)
    return p4est_quadrant_is_equal(q1,q2)
  elseif (Dc==3)
    return p8est_quadrant_is_equal(q1, q2)
  end
end

function _pXest_quadrant_is_parent(::Type{Val{Dc}}, q1, q2) where Dc
  if (Dc==2)
    return p4est_quadrant_is_parent(q1,q2)
  elseif (Dc==3)
    return p8est_quadrant_is_parent(q1,q2)
  end
end 

function _pXest_quadrant_compare(::Type{Val{Dc}}, q1, q2) where Dc
  if (Dc==2)
    return p4est_quadrant_compare(q1,q2)
  elseif (Dc==3)
    return p8est_quadrant_compare(q1,q2)
  end
end 

function _p4est_compute_migration_control_data(::Type{Val{Dc}},ptr_pXest_old,ptr_pXest_new) where Dc
  pXest_old   = ptr_pXest_old[]
  pXest_new   = ptr_pXest_new[]
  num_trees   = Cint(pXest_old.connectivity[].num_trees)
  my_rank     = pXest_old.mpirank
  ranks_count = Dict{Int,Int}()
  lst_ranks   = Int[]
  old2new     = Vector{Int}(undef,pXest_old.local_num_quadrants)
  current_old_quad_index = 1

  for itree = 0:num_trees-1
    tree = _pXest_tree_array_index(Val{Dc},pXest_old.trees,itree)[]
    num_quads = Cint(tree.quadrants.elem_count)

    for iquad = 0:num_quads-1
      q = _pXest_quadrant_array_index(Val{Dc},tree.quadrants, iquad)
      new_rank = _pXest_comm_find_owner(Val{Dc},ptr_pXest_new,itree,q,0)
      if (new_rank != my_rank)
        if (!(new_rank+1 in keys(ranks_count)))
          push!(lst_ranks,new_rank+1)
          ranks_count[new_rank+1] = 0
        end
        ranks_count[new_rank+1] += 1
        old2new[current_old_quad_index] = 0
      else
        current_new_quad_index = 1
        new_tree = _pXest_tree_array_index(Val{Dc},pXest_new.trees,pXest_new.first_local_tree)[]
        for t = pXest_new.first_local_tree:pXest_new.last_local_tree
          new_tree = _pXest_tree_array_index(Val{Dc},pXest_new.trees,t)[]
          if t == itree
            break
          end
          current_new_quad_index += Cint(new_tree.quadrants.elem_count)
        end
        found = false
        num_quads_new = Cint(new_tree.quadrants.elem_count)
        for iquad_new = 0:num_quads_new-1
          q_new = _pXest_quadrant_array_index(Val{Dc},new_tree.quadrants, iquad_new)
          found = _pXest_quadrant_is_equal(Val{Dc},q,q_new)!=0
          if found
            break
          end
          current_new_quad_index += 1
        end
        Gridap.Helpers.@check found
        old2new[current_old_quad_index] = current_new_quad_index
      end
      current_old_quad_index += 1
    end
  end

  local_ids    = [i for i=1:length(old2new) if old2new[i]==0]
  ptr_ranks    = Vector{Int32}(undef,length(ranks_count)+1)
  ptr_ranks[1] = 1
  for (i,rank) in enumerate(lst_ranks)
     ptr_ranks[i+1]=ptr_ranks[i]+ranks_count[rank]
  end

  lst_ranks,PartitionedArrays.JaggedArray(local_ids,ptr_ranks),old2new
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
function GridapDistributed.redistribute(model::OctreeDistributedDiscreteModel{Dc,Dp}, 
                                        parts_redistributed_model=model.parts) where {Dc,Dp}
  parts = (parts_redistributed_model === model.parts) ? model.parts : parts_redistributed_model
  comm  = parts.comm
  if (GridapDistributed.i_am_in(model.parts.comm) || GridapDistributed.i_am_in(parts.comm))
    if (parts_redistributed_model !== model.parts)
      A=is_included(model.parts,parts_redistributed_model)
      B=is_included(parts_redistributed_model,model.parts)
      @assert A || B
    end
    if (parts_redistributed_model===model.parts || A)
      _redistribute_parts_subseteq_parts_redistributed(model,parts_redistributed_model)
    else
      _redistribute_parts_supset_parts_redistributed(model, parts_redistributed_model)
    end
  else
    VoidOctreeDistributedDiscreteModel(model,model.parts), nothing
  end
end

function _redistribute_parts_subseteq_parts_redistributed(model::OctreeDistributedDiscreteModel{Dc,Dp}, parts_redistributed_model) where {Dc,Dp}
  parts = (parts_redistributed_model === model.parts) ? model.parts : parts_redistributed_model
  if (parts_redistributed_model === model.parts)
    ptr_pXest_old = model.ptr_pXest
  else
    ptr_pXest_old = _p4est_to_new_comm(model.ptr_pXest,
                                       model.ptr_pXest_connectivity,
                                       model.parts.comm,
                                       parts.comm)
  end
  ptr_pXest_new = pXest_copy(Val{Dc}, ptr_pXest_old)
  pXest_partition!(Val{Dc}, ptr_pXest_new)

  # Compute RedistributeGlue
  parts_snd, lids_snd, old2new = _p4est_compute_migration_control_data(Val{Dc},ptr_pXest_old,ptr_pXest_new)
  parts_rcv, lids_rcv, new2old = _p4est_compute_migration_control_data(Val{Dc},ptr_pXest_new,ptr_pXest_old)

  lids_rcv, parts_rcv, lids_snd, parts_snd, old2new, new2old =
       _to_pdata(parts, lids_rcv, parts_rcv, lids_snd, parts_snd, old2new, new2old)

  glue = GridapDistributed.RedistributeGlue(parts,model.parts,parts_rcv,parts_snd,lids_rcv,lids_snd,old2new,new2old)

  # Extract ghost and lnodes
  ptr_pXest_ghost  = setup_pXest_ghost(Val{Dc}, ptr_pXest_new)
  ptr_pXest_lnodes = setup_pXest_lnodes_nonconforming(Val{Dc}, ptr_pXest_new, ptr_pXest_ghost)

  # Build fine-grid mesh
  fmodel, non_conforming_glue = setup_non_conforming_distributed_discrete_model(Val{Dc},
                                            parts,
                                            model.coarse_model,
                                            model.ptr_pXest_connectivity,
                                            ptr_pXest_new,
                                            ptr_pXest_ghost,
                                            ptr_pXest_lnodes)

  pXest_lnodes_destroy(Val{Dc},ptr_pXest_lnodes)
  pXest_ghost_destroy(Val{Dc},ptr_pXest_ghost)

  red_model = OctreeDistributedDiscreteModel(Dc,Dp,
                                             parts,
                                             fmodel,
                                             non_conforming_glue,
                                             model.coarse_model,
                                             model.ptr_pXest_connectivity,
                                             ptr_pXest_new,
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
    first_global_quadrant=Int64((Float64(N)*Float64(psub-1))/(Float64(Psub)))
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

  ptr_pXest_old=pXest_copy(Val{Dc},model.ptr_pXest)
  pXest_partition_given!(Val{Dc}, ptr_pXest_old, num_cells_per_part)

  # p4est_vtk_write_file(ptr_pXest_old, C_NULL, "ptr_pXest_old")

  # ptr_pXest_old is distributed over supset_comm
  # once created, ptr_pXest_new is distributed over subset_comm
  ptr_pXest_new = _p4est_to_new_comm(ptr_pXest_old,
                                     model.ptr_pXest_connectivity,
                                     supset_comm,
                                     subset_comm)

  # Compute RedistributeGlue
  parts_snd, lids_snd, old2new =
      _p4est_compute_migration_control_data(Val{Dc},model.ptr_pXest,ptr_pXest_old)
  parts_rcv, lids_rcv, new2old =
      _p4est_compute_migration_control_data(Val{Dc},ptr_pXest_old,model.ptr_pXest)

  pXest_destroy(Val{Dc},ptr_pXest_old)

  lids_rcv, parts_rcv, lids_snd, parts_snd, old2new, new2old =
       _to_pdata(model.parts, lids_rcv, parts_rcv, lids_snd, parts_snd, old2new, new2old)

  glue = GridapDistributed.RedistributeGlue(parts_redistributed_model,model.parts,
                                            parts_rcv,parts_snd,lids_rcv,lids_snd,
                                            old2new,new2old)

  if (GridapDistributed.i_am_in(subset_comm))
    # p4est_vtk_write_file(ptr_pXest_new, C_NULL, "ptr_pXest_new")

    # Extract ghost and lnodes
    ptr_pXest_ghost  = setup_pXest_ghost(Val{Dc}, ptr_pXest_new)
    ptr_pXest_lnodes = setup_pXest_lnodes(Val{Dc}, ptr_pXest_new, ptr_pXest_ghost)

    # # Build fine-grid mesh
    fmodel = setup_distributed_discrete_model(Val{Dc},
                                              parts_redistributed_model,
                                              model.coarse_model,
                                              model.ptr_pXest_connectivity,
                                              ptr_pXest_new,
                                              ptr_pXest_ghost,
                                              ptr_pXest_lnodes)

    pXest_lnodes_destroy(Val{Dc},ptr_pXest_lnodes)
    pXest_ghost_destroy(Val{Dc},ptr_pXest_ghost)

    non_conforming_glue = _create_conforming_model_non_conforming_glue(fmodel)

    red_model = OctreeDistributedDiscreteModel(Dc,Dp,
                                               parts_redistributed_model,
                                               fmodel,
                                               non_conforming_glue,
                                               model.coarse_model,
                                               model.ptr_pXest_connectivity,
                                               ptr_pXest_new,
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
                           
const p8est_subface_to_hanging_edges_within_subface = 
[ 
  1 3;
  1 2;
  0 3;
  0 2;
]   

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


# To add to P4est_wrapper.jl library
# I just translated this function to Julia from its p4est counterpart
# We cannot call it directly because it is declared as static within p4est,
# and thus it does not belong to the ABI of the dynamic library object.

# /** Decode the face_code into hanging face information.
#  *
#  * This is mostly for demonstration purposes.  Applications probably will
#  * integrate it into their own loop over the face for performance reasons.
#  *
#  * \param[in] face_code as in the p4est_lnodes_t structure.
#  * \param[out] hanging face: if there are hanging faces,
#  *             hanging_face = -1 if the face is not hanging,
#  *                          = 0 if the face is the first half,
#  *                          = 1 if the face is the second half.
#  *             note: not touched if there are no hanging faces.
#  * \return              true if any face is hanging, false otherwise.
#  */

const p4est_corner_faces = [0 2; 1 2; 0 3; 1 3]
const p4est_corner_face_corners = [0 -1 0 -1; -1 0 1 -1; 1 -1 -1 0; -1 1 -1 1]
function p4est_lnodes_decode(face_code, hanging_face)
  @assert face_code >= 0
  if (face_code != 0)
    c = face_code & 0x03
    work = face_code >> 2
    hanging_face .= -1
    for i = 0:1
      f = p4est_corner_faces[c+1, i+1]
      hanging_face[f+1] = (work & 0x01) != 0 ? p4est_corner_face_corners[c+1, f+1] : -1
      work >>= 1
    end
    return 1
  else
    return 0
  end
end

const p8est_corner_faces = [0 2 4; 1 2 4; 0 3 4; 1 3 4; 0 2 5; 1 2 5; 0 3 5; 1 3 5]

const p8est_face_edges = [ 4 6 8 10; 5 7 9 11; 0 2 8 9; 1 3 10 11; 0 1 4 5; 2 3 6 7]

const p8est_corner_face_corners = [0 -1  0 -1  0 -1; -1  0  1 -1  1 -1 ; 1 -1 -1  0  2 -1 ; -1  1 -1  1  3 -1 ;
                                   2 -1  2 -1 -1  0; -1  2  3 -1 -1  1 ; 3 -1 -1  2 -1  2 ; -1  3 -1  3 -1  3 ]

const p8est_corner_edges = [ 0 4 8; 0 5 9; 1 4 10; 1 5 11; 2 6 8; 2 7 9; 3 6 10; 3 7 11 ]

# To add to p8est_wrapper.jl library
# I just translated this function to Julia from its p4est counterpart
# We cannot call it directly because it is declared as static within p4est,
# and thus it does not belong to the ABI of the dynamic library object.
function p8est_lnodes_decode(face_code,
                             hanging_face,
                             hanging_edge)
  @assert face_code >= 0
  if (face_code!=0)
    c = face_code & 0x0007
    work = face_code >> 3
    hanging_face .= -1
    hanging_edge .= -1
    cwork = c
    for i=0:2
      if ((work & 0x0001)!=0)
        f = p8est_corner_faces[c+1,i+1]
        hanging_face[f+1] = p8est_corner_face_corners[c+1,f+1]
        for j=0:3
          e = p8est_face_edges[f+1,j+1]
          hanging_edge[e+1] = 4
        end
      end
      work >>= 1
    end
    for i=0:3
      if ((work & 0x0001)!=0)
        e = p8est_corner_edges[c+1,i+1]
        hanging_edge[e+1] = (hanging_edge[e+1] == -1) ? 0 : 2
        hanging_edge[e+1] += (cwork & 0x0001)
      end
      cwork >>= 1
      work  >>= 1
    end
    return 1
  else
    return 0
  end
end

function setup_non_conforming_distributed_discrete_model(::Type{Val{Dc}},
                                                         parts,
                                                         coarse_discrete_model,
                                                         ptr_pXest_connectivity,
                                                         ptr_pXest,
                                                         ptr_pXest_ghost,
                                                         ptr_pXest_lnodes) where Dc

  cell_prange = setup_cell_prange(Val{Dc}, parts, ptr_pXest, ptr_pXest_ghost)

  gridap_cell_faces,
  non_conforming_glue=
    generate_cell_faces_and_non_conforming_glue(Val{Dc},ptr_pXest_lnodes, cell_prange)


  nlvertices = map(non_conforming_glue) do ncglue
    ncglue.num_regular_faces[1]+ncglue.num_hanging_faces[1]
  end

  node_coordinates=generate_node_coordinates(Val{Dc},
                                             gridap_cell_faces[1],
                                             nlvertices,
                                             ptr_pXest_connectivity,
                                             ptr_pXest,
                                             ptr_pXest_ghost)

  grid,topology=generate_grid_and_topology(Val{Dc},
                                           gridap_cell_faces[1],
                                           nlvertices,
                                           node_coordinates)

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

  face_labeling=generate_face_labeling(parts,
                                       cell_prange,
                                       coarse_discrete_model,
                                       grid,
                                       topology,
                                       ptr_pXest,
                                       ptr_pXest_ghost)

  _set_hanging_labels!(face_labeling,non_conforming_glue)

  discretemodel=map(grid,topology,face_labeling) do grid, topology, face_labeling
    Gridap.Geometry.UnstructuredDiscreteModel(grid,topology,face_labeling)
  end
  GridapDistributed.DistributedDiscreteModel(discretemodel,cell_prange), non_conforming_glue
end

function _set_hanging_labels!(face_labeling,non_conforming_glue)
  max_entity_ids = map(face_labeling) do face_labeling
    max_entity_id = typemin(eltype(first(face_labeling.d_to_dface_to_entity))) 
    for i=1:length(face_labeling.d_to_dface_to_entity)
      max_entity_id=max(maximum(face_labeling.d_to_dface_to_entity[i]),max_entity_id)
    end
    max_entity_id
  end
  max_entity_id = reduction(max,
                            max_entity_ids,
                            destination=:all,
                            init=zero(eltype(max_entity_ids)))
  
  hanging_entitity_ids = Dict{Int,Bool}()
  map(max_entity_id,
      face_labeling,
      non_conforming_glue) do max_entity_id,face_labeling,ncglue 
    for i=1:length(ncglue.num_hanging_faces)
      num_regular_faces_i = ncglue.num_regular_faces[i]
      num_hanging_faces_i = ncglue.num_hanging_faces[i]
      for j=num_regular_faces_i+1:num_regular_faces_i+num_hanging_faces_i
        hanging_entity_id = max_entity_id + face_labeling.d_to_dface_to_entity[i][j]
        face_labeling.d_to_dface_to_entity[i][j]=hanging_entity_id
        hanging_entitity_ids[hanging_entity_id]=true
      end      
    end
  end
  map(face_labeling) do face_labeling 
    add_tag!(face_labeling,"hanging",collect(keys(hanging_entitity_ids)))
  end
end 

function _build_map_from_faces_to_cell_lface(vnodes, element_nodes, face_code)
  n_cell_faces    = num_cell_faces(Val{2})
  hanging_face = Vector{Cint}(undef, n_cell_faces)

  # Build a map from faces to (cell,lface)
  p4est_gface_to_gcell_p4est_lface = Dict{Int,Tuple{Int,Int}}()
  for cell = 1:length(face_code)
    start = (cell - 1) * vnodes + 1
    p4est_cell_faces = view(element_nodes, start:start+n_cell_faces-1)
    has_hanging = p4est_lnodes_decode(face_code[cell], hanging_face)
    if (has_hanging==0)
      for (lface, gface) in enumerate(p4est_cell_faces)
        p4est_gface_to_gcell_p4est_lface[gface] = (cell, lface)
      end 
    else
      for (lface, half) in enumerate(hanging_face)
        # Current face is NOT hanging
        if (half == -1)
          gface = p4est_cell_faces[lface]
          p4est_gface_to_gcell_p4est_lface[gface] = (cell, lface)
        end   
      end
    end
  end
  p4est_gface_to_gcell_p4est_lface
end 

function _build_map_from_faces_edges_to_cell_lface_ledge(vnodes, element_nodes, face_code)
    n_cell_faces    = num_cell_faces(Val{3})
    n_cell_edges    = num_cell_edges(Val{3})

    hanging_face = Vector{Cint}(undef, n_cell_faces)
    hanging_edge = Vector{Cint}(undef, n_cell_edges)

    # Build a map from faces to (cell,lface)
    p4est_gface_to_gcell_p4est_lface = Dict{Int,Tuple{Int,Int}}()
    p4est_gedge_to_gcell_p4est_ledge = Dict{Int,Tuple{Int,Int}}()
    for cell = 1:length(face_code)
      start = (cell - 1) * vnodes + 1
      p4est_cell_faces = view(element_nodes, start:start+n_cell_faces-1)
      p4est_cell_edges = view(element_nodes, start+n_cell_faces:start+n_cell_faces+n_cell_edges-1)

      
      has_hanging = p8est_lnodes_decode(face_code[cell], hanging_face, hanging_edge)
      if (has_hanging==0)
        for (lface, gface) in enumerate(p4est_cell_faces)
          p4est_gface_to_gcell_p4est_lface[gface] = (cell, lface)
        end 
        for (ledge, gedge) in enumerate(p4est_cell_edges)
          p4est_gedge_to_gcell_p4est_ledge[gedge] = (cell, ledge)
        end 
      else
        for (lface, half) in enumerate(hanging_face)
          # Current face is NOT hanging
          if (half == -1)
            gface = p4est_cell_faces[lface]
            p4est_gface_to_gcell_p4est_lface[gface] = (cell, lface)
          end   
        end
        for (ledge, half) in enumerate(hanging_edge)
          # Current edge is NOT hanging
          if (half == -1)
            gedge = p4est_cell_edges[ledge]
            p4est_gedge_to_gcell_p4est_ledge[gedge] = (cell, ledge)
          end   
        end
      end
    end
    p4est_gface_to_gcell_p4est_lface, p4est_gedge_to_gcell_p4est_ledge
end

function pXest_2_gridap_vertex(::Type{Val{Dc}}) where Dc
  Gridap.Arrays.IdentityVector(num_cell_vertices(Val{Dc}))
end

function p8est_2_gridap_edge()
  Gridap.Arrays.IdentityVector(num_cell_edges(Val{3}))
end

function pXest_2_gridap_facet(::Type{Val{Dc}}) where Dc
  if (Dc==2)
    GridapP4est.P4EST_2_GRIDAP_FACET_2D
  else
    @assert Dc==3 
    GridapP4est.P4EST_2_GRIDAP_FACET_3D
  end 
end

function hanging_lvertex_within_face_2d(half)
  half == 0 ? 1 : 0
end

function hanging_lvertex_within_face_3d(half)
  if (half==0)
    return 3
  elseif (half==1)
    return 2 
  elseif (half==2)
    return 1 
  elseif (half==3)
    return 0
  end 
end

function hanging_lvertex_within_edge(half)
  if (half==0 || half==2)
    return 1 
  elseif (half==1 || half==3)
    return 0 
  end 
  @assert false
end

function regular_lvertex_within_face(half)
  return half
end

function regular_lvertex_within_edge(half)
  if (half==0 || half==2)
    return 0 
  elseif (half==1 || half==3)
    return 1 
  end 
  @assert false
end


function generate_cell_faces_and_non_conforming_glue(::Type{Val{Dc}}, 
                                                     ptr_pXest_lnodes, 
                                                     cell_prange) where Dc
  
  n_cell_vertices = num_cell_vertices(Val{Dc})
  n_cell_edges    = num_cell_edges(Val{Dc})
  n_cell_faces    = num_cell_faces(Val{Dc})
  
  lnodes = ptr_pXest_lnodes[]
  element_nodes = unsafe_wrap(Array, lnodes.element_nodes, lnodes.vnodes * lnodes.num_local_elements)
  face_code = unsafe_wrap(Array, lnodes.face_code, lnodes.num_local_elements)
  hanging_face = Vector{Cint}(undef, n_cell_faces)
  face_code_with_ghosts = map(partition(cell_prange)) do indices
      @assert length(face_code)==own_length(indices)
      @assert own_length(indices)==lnodes.num_local_elements
      face_code_with_ghosts=similar(face_code, local_length(indices))
      face_code_with_ghosts[1:own_length(indices)] .= face_code
      face_code_with_ghosts
  end

  cache_face_code=fetch_vector_ghost_values_cache(face_code_with_ghosts, partition(cell_prange))
  fetch_vector_ghost_values!(face_code_with_ghosts, cache_face_code) |> wait

  element_nodes_with_ghosts = map(partition(cell_prange)) do indices
    nonlocal_nodes = unsafe_wrap(Array, lnodes.nonlocal_nodes, lnodes.num_local_nodes-lnodes.owned_count)
    element_nodes_with_ghosts_data=similar(element_nodes, local_length(indices)*lnodes.vnodes)
    for (i,node) in enumerate(element_nodes)
      if (node<lnodes.owned_count)
        element_nodes_with_ghosts_data[i] = lnodes.global_offset+node
      else
        element_nodes_with_ghosts_data[i] = nonlocal_nodes[node-lnodes.owned_count+1]
      end 
    end     
    element_nodes_with_ghosts_ptrs = [i for i=1:lnodes.vnodes:length(element_nodes_with_ghosts_data)+1]
    PArrays.JaggedArray(element_nodes_with_ghosts_data,element_nodes_with_ghosts_ptrs)
  end
  cache_element_nodes_with_ghosts=fetch_vector_ghost_values_cache(element_nodes_with_ghosts, partition(cell_prange))
  fetch_vector_ghost_values!(element_nodes_with_ghosts, cache_element_nodes_with_ghosts) |> wait

  map(element_nodes_with_ghosts,face_code_with_ghosts,partition(cell_prange)) do element_nodes_with_ghosts, face_code_with_ghosts, indices
    @debug "ENDES[$(part_id(indices))]: $(element_nodes_with_ghosts.data)"
    @debug "FCODS[$(part_id(indices))]: $(face_code_with_ghosts)"
  end

  if (Dc==2)
    hanging_lvertex_within_face=hanging_lvertex_within_face_2d
    pXest_face_corners = p4est_face_corners 
  else
    hanging_lvertex_within_face=hanging_lvertex_within_face_3d
    pXest_face_corners = p8est_face_corners
  end 

  if (Dc==3)
    hanging_edge = Vector{Cint}(undef, n_cell_edges)
  end

  num_regular_faces,
  num_hanging_faces,
  gridap_cell_faces,
  hanging_faces_glue = 
      map(partition(cell_prange),
          element_nodes_with_ghosts,
          face_code_with_ghosts) do indices, 
                                    element_nodes_with_ghosts, 
                                    face_code_with_ghosts
    @assert local_length(indices)==length(face_code_with_ghosts) 
 
    num_local_elements = local_length(indices)
    num_regular_faces = Vector{Int}(undef, Dc)
    num_hanging_faces = Vector{Int}(undef, Dc)

    # Count regular vertices
    num_regular_faces[1] = 0
    regular_vertices_p4est_to_gridap = Dict{Int,Int}()

    num_regular_faces[Dc] = 0
    regular_faces_p4est_to_gridap = Dict{Int,Int}()

    if (Dc==2)
      p4est_gface_to_gcell_p4est_lface = 
         _build_map_from_faces_to_cell_lface(lnodes.vnodes, element_nodes_with_ghosts.data, face_code_with_ghosts)
    else 
      p4est_gface_to_gcell_p4est_lface, 
         p4est_gedge_to_gcell_p4est_ledge = 
           _build_map_from_faces_edges_to_cell_lface_ledge(lnodes.vnodes, 
                                                            element_nodes_with_ghosts.data, 
                                                            face_code_with_ghosts)
    end 

    PXEST_2_GRIDAP_VERTEX = pXest_2_gridap_vertex(Val{Dc})
    PXEST_2_GRIDAP_FACE   = pXest_2_gridap_facet(Val{Dc})
    PXEST_2_GRIDAP_EDGE   = p8est_2_gridap_edge()

    n = local_length(indices)
    gridap_cell_vertices_ptrs = Vector{Int32}(undef,n+1)
    gridap_cell_faces_ptrs = Vector{Int32}(undef,n+1)
    gridap_cell_vertices_ptrs[1]=1
    gridap_cell_faces_ptrs[1]=1

    hanging_vertices_pairs_to_owner_face = Dict{Tuple{Int,Int},Int}()
    hanging_faces_pairs_to_owner_face = Dict{Tuple{Int,Int},Tuple{Int,Int}}()

    for i=1:n
      gridap_cell_vertices_ptrs[i+1]=gridap_cell_vertices_ptrs[i]+n_cell_vertices
      gridap_cell_faces_ptrs[i+1]=gridap_cell_faces_ptrs[i]+n_cell_faces
    end

    gridap_cell_vertices_data = Vector{Int}(undef, num_local_elements * n_cell_vertices)
    gridap_cell_vertices_data .= -1

    gridap_cell_faces_data = Vector{Int}(undef, num_local_elements * n_cell_faces)
    gridap_cell_faces_data .= -1

    if (Dc==3)
      num_regular_faces[2] = 0

      gridap_cell_edges_ptrs = Vector{Int32}(undef,n+1)
      gridap_cell_edges_ptrs[1]=1
      for i=1:n
        gridap_cell_edges_ptrs[i+1]=gridap_cell_edges_ptrs[i]+n_cell_edges
      end
      gridap_cell_edges_data = Vector{Int}(undef, num_local_elements * n_cell_edges)
      gridap_cell_edges_data .= -1
      hanging_edges_cell_ledge_to_owner_face_half = Dict{Tuple{Int,Int},Tuple{Int,Int}}()
      owner_edge_subedge_to_cell_ledge = Dict{Tuple{Int,Int},Tuple{Int,Int}}()
      hanging_vertices_pairs_to_owner_edge = Dict{Tuple{Int,Int},Int}()
      regular_edges_p4est_to_gridap = Dict{Int,Int}()
    end 

    for cell = 1:num_local_elements
      start                 = (cell - 1) * lnodes.vnodes + 1
      start_gridap_vertices = (cell - 1) * n_cell_vertices
      start_gridap_faces    = (cell - 1) * n_cell_faces

      p4est_cell_faces = view(element_nodes_with_ghosts.data, start:start+n_cell_faces-1)
      p4est_cell_vertices = view(element_nodes_with_ghosts.data, 
                                 start+n_cell_faces+n_cell_edges:start+n_cell_faces+n_cell_edges+n_cell_vertices-1)

      gridap_cell_vertices = view(gridap_cell_vertices_data,
        start_gridap_vertices+1:start_gridap_vertices+n_cell_vertices)
      gridap_cell_faces = view(gridap_cell_faces_data,
        start_gridap_faces+1:start_gridap_faces+n_cell_faces)

      if (Dc==2)  
        has_hanging = p4est_lnodes_decode(face_code_with_ghosts[cell], hanging_face)
      else
        has_hanging = p8est_lnodes_decode(face_code_with_ghosts[cell], hanging_face, hanging_edge)
        start_gridap_edges = (cell-1)*n_cell_edges
        gridap_cell_edges = view(gridap_cell_edges_data, start_gridap_edges+1:start_gridap_edges+n_cell_edges)
        p4est_cell_edges = view(element_nodes_with_ghosts.data, 
                                start+n_cell_faces:start+n_cell_faces+n_cell_edges-1)
      end
      if has_hanging == 0
        # All vertices/edges/faces of the current cell are regular 
        # Process vertices
        for (p4est_lvertex, p4est_gvertex) in enumerate(p4est_cell_vertices)
          num_regular_faces[1] =
            process_current_face!(gridap_cell_vertices,
              regular_vertices_p4est_to_gridap,
              num_regular_faces[1],
              p4est_cell_vertices,
              p4est_lvertex,
              p4est_gvertex,
              PXEST_2_GRIDAP_VERTEX)
        end
        
        if (Dc==3)
          for (p4est_ledge, p4est_gedge) in enumerate(p4est_cell_edges)
            num_regular_faces[2] =
              process_current_face!(gridap_cell_edges,
                regular_edges_p4est_to_gridap,
                num_regular_faces[2],
                p4est_cell_edges,
                p4est_ledge,
                p4est_gedge,
                PXEST_2_GRIDAP_EDGE)
          end
        end 

        # Process faces
        for (p4est_lface, p4est_gface) in enumerate(p4est_cell_faces)
          num_regular_faces[Dc] =
            process_current_face!(gridap_cell_faces,
              regular_faces_p4est_to_gridap,
              num_regular_faces[Dc],
              p4est_cell_faces,
              p4est_lface,
              p4est_gface,
              PXEST_2_GRIDAP_FACE)
        end
      else
        # "Touch" hanging vertices before processing current cell
        # This is required as we dont have any means to detect 
        # a hanging vertex from a non-hanging face
        for (p4est_lface, half) in enumerate(hanging_face)
          if (half != -1)
            hanging_vertex_lvertex_within_face = hanging_lvertex_within_face(half)
            p4est_lvertex = pXest_face_corners[p4est_lface,
                                               hanging_vertex_lvertex_within_face+1]
            gridap_cell_vertices[PXEST_2_GRIDAP_VERTEX[p4est_lvertex+1]] = hanging_vertex_code
          end 
        end

        if (Dc==3)
          for (p4est_ledge, half) in enumerate(hanging_edge)
            if (half != -1 && half !=4)
              hanging_vertex_lvertex_within_edge = hanging_lvertex_within_edge(half)
              p4est_lvertex = p8est_edge_corners[p4est_ledge,
                                                 hanging_vertex_lvertex_within_edge+1]
              gridap_cell_vertices[PXEST_2_GRIDAP_VERTEX[p4est_lvertex+1]] = hanging_vertex_code
            end 
          end
        end 

        # Current cell has at least one hanging face 
        for (p4est_lface, half) in enumerate(hanging_face)
          # Current face is NOT hanging
          if (half == -1)
            # Process vertices on the boundary of p4est_lface
            for p4est_lvertex in pXest_face_corners[p4est_lface, :]
              p4est_gvertex = p4est_cell_vertices[p4est_lvertex+1]
              if (gridap_cell_vertices[p4est_lvertex+1] != hanging_vertex_code)
                num_regular_faces[1] =
                  process_current_face!(gridap_cell_vertices,
                    regular_vertices_p4est_to_gridap,
                    num_regular_faces[1],
                    p4est_cell_vertices,
                    p4est_lvertex + 1,
                    p4est_gvertex,
                    PXEST_2_GRIDAP_VERTEX)
              end
            end
            # Process non-hanging face
            p4est_gface = p4est_cell_faces[p4est_lface]
            num_regular_faces[Dc] =
              process_current_face!(gridap_cell_faces,
                regular_faces_p4est_to_gridap,
                num_regular_faces[Dc],
                p4est_cell_faces,
                p4est_lface,
                p4est_gface,
                PXEST_2_GRIDAP_FACE)
          else # Current face is hanging
            # Identify regular vertex and hanging vertex 
            # Repeat code above for regular vertex 
            # Special treatment for hanging vertex 
            regular_vertex_lvertex_within_face = regular_lvertex_within_face(half)
            hanging_vertex_lvertex_within_face = hanging_lvertex_within_face(half)

            # Process regular vertex
            p4est_regular_lvertex = pXest_face_corners[p4est_lface, regular_vertex_lvertex_within_face+1]
            p4est_gvertex = p4est_cell_vertices[p4est_regular_lvertex+1]
            num_regular_faces[1] =
              process_current_face!(gridap_cell_vertices,
                regular_vertices_p4est_to_gridap,
                num_regular_faces[1],
                p4est_cell_vertices,
                p4est_regular_lvertex + 1,
                p4est_gvertex,
                PXEST_2_GRIDAP_VERTEX)
            
            # Process hanging vertex
            p4est_hanging_lvertex = pXest_face_corners[p4est_lface, hanging_vertex_lvertex_within_face+1]
            owner_face = p4est_cell_faces[p4est_lface]
            hanging_vertices_pairs_to_owner_face[(cell, PXEST_2_GRIDAP_VERTEX[p4est_hanging_lvertex+1])] = owner_face
            
            # Process hanging face
            hanging_faces_pairs_to_owner_face[(cell, PXEST_2_GRIDAP_FACE[p4est_lface])] = (owner_face,half+1)

            if (Dc==3)
              for (i,ledge_within_face) in enumerate(p8est_subface_to_hanging_edges_within_subface[half+1,:])
                p4est_ledge=p8est_face_edges[p4est_lface,ledge_within_face+1]
                gridap_ledge = PXEST_2_GRIDAP_EDGE[p4est_ledge+1]
                # Identify the two edges which are hanging within the face
                hanging_edges_cell_ledge_to_owner_face_half[(cell, gridap_ledge)] =
                    (owner_face,-p8est_subface_to_hanging_edges_within_face[half+1,i])
                gridap_cell_edges[gridap_ledge] = hanging_edge_from_face_code
              end 
            end 

          end
        end


        if (Dc==3)
          for (p4est_ledge, half) in enumerate(hanging_edge)
            # Current edge is NOT hanging
            if (half == -1)
              # Process vertices on the boundary of p4est_ledge
              for p4est_lvertex in p8est_edge_corners[p4est_ledge, :]
                p4est_gvertex = p4est_cell_vertices[p4est_lvertex+1]
                if (gridap_cell_vertices[p4est_lvertex+1] != hanging_vertex_code)
                  num_regular_faces[1] =
                    process_current_face!(gridap_cell_vertices,
                      regular_vertices_p4est_to_gridap,
                      num_regular_faces[1],
                      p4est_cell_vertices,
                      p4est_lvertex + 1,
                      p4est_gvertex,
                      PXEST_2_GRIDAP_VERTEX)
                end
              end
              # Process non-hanging edge
              p4est_gedge = p4est_cell_edges[p4est_ledge]
              num_regular_faces[2] =
                process_current_face!(gridap_cell_edges,
                  regular_edges_p4est_to_gridap,
                  num_regular_faces[2],
                  p4est_cell_edges,
                  p4est_ledge,
                  p4est_gedge,
                  PXEST_2_GRIDAP_EDGE)
            else # Current edge is hanging
              if ( gridap_cell_edges[PXEST_2_GRIDAP_EDGE[p4est_ledge]] != hanging_edge_from_face_code )
                # The present hanging edge cannot be within a coarser face
                @assert half != 4 

                # # Identify regular vertex and hanging vertex 
                # # Repeat code above for regular vertex 
                # # Special treatment for hanging vertex 
                regular_vertex_lvertex_within_edge = regular_lvertex_within_edge(half)
                hanging_vertex_lvertex_within_edge = hanging_lvertex_within_edge(half)

                # # Process regular vertex
                p4est_regular_lvertex = p8est_edge_corners[p4est_ledge, regular_vertex_lvertex_within_edge+1]
                p4est_gvertex = p4est_cell_vertices[p4est_regular_lvertex+1]
                
                num_regular_faces[1] =
                  process_current_face!(gridap_cell_vertices,
                    regular_vertices_p4est_to_gridap,
                    num_regular_faces[1],
                    p4est_cell_vertices,
                    p4est_regular_lvertex + 1,
                    p4est_gvertex,
                    PXEST_2_GRIDAP_VERTEX)
                
                # Process hanging vertex
                p4est_hanging_lvertex = p8est_edge_corners[p4est_ledge, hanging_vertex_lvertex_within_edge+1]
                p4est_owner_edge = p4est_cell_edges[p4est_ledge]
                hanging_vertices_pairs_to_owner_edge[(cell, 
                                                      PXEST_2_GRIDAP_VERTEX[p4est_hanging_lvertex+1])] = p4est_owner_edge

                # Process hanging edge
                subedge = regular_vertex_lvertex_within_edge+1
                owner_edge_subedge_pair=(p4est_owner_edge,subedge)
                gridap_ledge=PXEST_2_GRIDAP_EDGE[p4est_ledge]
                hanging_edges_cell_ledge_to_owner_face_half[(cell, gridap_ledge)] = owner_edge_subedge_pair
                if (!haskey(owner_edge_subedge_to_cell_ledge,owner_edge_subedge_pair))
                  owner_edge_subedge_to_cell_ledge[owner_edge_subedge_pair] = (cell,gridap_ledge)
                end
              end 
            end 
        end
      end
    end
  end 

    function is_ghost(cell)
      cell>own_length(indices)
    end

    # Go over all touched hanging faces and start 
    # assigning IDs from the last num_regular_faces ID
    # For each hanging face, keep track of (owner_cell,lface)
    # Go over all hanging faces 
    # Detect if the owner face is in a ghost cell. 
    # If not in a ghost cell or touched 
    # Else 
    #   The current face becomes a regular face 
    # end 
    hanging_faces_owner_cell_and_lface =
      Vector{Tuple{Int,Int,Int}}(undef, length(keys(hanging_faces_pairs_to_owner_face)))
    num_hanging_faces[Dc] = 0
    for key in keys(hanging_faces_pairs_to_owner_face)
      (cell, lface) = key
      (owner_p4est_gface, half) = hanging_faces_pairs_to_owner_face[key]
      num_hanging_faces[Dc] += 1
      start_gridap_faces = (cell - 1) * n_cell_faces
      gridap_cell_faces_data[start_gridap_faces+lface] = num_regular_faces[Dc] + num_hanging_faces[Dc]
      if (!(is_ghost(cell)) || haskey(regular_faces_p4est_to_gridap,owner_p4est_gface))  
        @assert haskey(regular_faces_p4est_to_gridap,owner_p4est_gface)
        (owner_cell, p4est_lface) = p4est_gface_to_gcell_p4est_lface[owner_p4est_gface]
        hanging_faces_owner_cell_and_lface[num_hanging_faces[Dc]] =
          (owner_cell, n_cell_vertices+n_cell_edges+PXEST_2_GRIDAP_FACE[p4est_lface], half)
      else
        # Glue info cannot be computed for this hanging face
        hanging_faces_owner_cell_and_lface[num_hanging_faces[Dc]] = (-1,-1,-1)
      end
    end

    @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))]  gridap_cell_faces_data: $(gridap_cell_faces_data)"


    # Go over all touched hanging vertices and start 
    # assigning IDs from the last num_regular_vertices ID
    # For each hanging vertex, keep track of (owner_cell,lface)
    num_hanging_faces[1] = 0
    hanging_vertices_owner_cell_and_lface = Tuple{Int,Int,Int}[]
    half=1
    owner_p4est_gface_to_hanging_vertex = Dict{Int,Int}()
    for key in keys(hanging_vertices_pairs_to_owner_face)
      (cell, lvertex) = key
      owner_p4est_gface = hanging_vertices_pairs_to_owner_face[key]
      if !(haskey(owner_p4est_gface_to_hanging_vertex, owner_p4est_gface))
        num_hanging_faces[1] += 1
        owner_p4est_gface_to_hanging_vertex[owner_p4est_gface] = num_hanging_faces[1]
        if (!is_ghost(cell) || (haskey(regular_faces_p4est_to_gridap,owner_p4est_gface)))
          (owner_cell, p4est_lface) = p4est_gface_to_gcell_p4est_lface[owner_p4est_gface]
          push!(hanging_vertices_owner_cell_and_lface,
            (owner_cell, n_cell_vertices+n_cell_edges+PXEST_2_GRIDAP_FACE[p4est_lface],half))
        else
          push!(hanging_vertices_owner_cell_and_lface,(-1, -1,-1))
        end
      end 
      start_gridap_vertices = (cell - 1) * n_cell_vertices
      gridap_cell_vertices_data[start_gridap_vertices+lvertex] = num_regular_faces[1] +
                                                          owner_p4est_gface_to_hanging_vertex[owner_p4est_gface] 
    end

    @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))]  gridap_cell_vertices_data: $(gridap_cell_vertices_data)"
      
    if (Dc==3)
      half=1
      owner_p4est_gedge_to_hanging_vertex = Dict{Int,Int}()
      for key in keys(hanging_vertices_pairs_to_owner_edge)
        (cell, lvertex) = key
        owner_p4est_gedge = hanging_vertices_pairs_to_owner_edge[key]
        if !(haskey(owner_p4est_gedge_to_hanging_vertex, owner_p4est_gedge))
          num_hanging_faces[1] += 1
          owner_p4est_gedge_to_hanging_vertex[owner_p4est_gedge] = num_hanging_faces[1]
          if (!is_ghost(cell) || (haskey(regular_edges_p4est_to_gridap,owner_p4est_gedge)))
              (owner_cell, p4est_ledge) = p4est_gedge_to_gcell_p4est_ledge[owner_p4est_gedge]
              push!(hanging_vertices_owner_cell_and_lface,
                      (owner_cell, n_cell_vertices+PXEST_2_GRIDAP_EDGE[p4est_ledge],half))
          else 
            push!(hanging_vertices_owner_cell_and_lface,(-1, -1,-1))
          end
        end
        start_gridap_vertices = (cell - 1) * n_cell_vertices
        gridap_cell_vertices_data[start_gridap_vertices+lvertex] = num_regular_faces[1] +
                                                owner_p4est_gedge_to_hanging_vertex[owner_p4est_gedge]
      end 


      # Go over all touched hanging edges and start 
      # assigning IDs from the last num_regular_edge ID
      # For each hanging edge, keep track of (owner_cell,lface/ledge)
      hanging_edges_owner_cell_and_lface = Tuple{Int,Int,Int}[]
      owner_p4est_gface_half_to_hanging_edge = Dict{Tuple{Int,Int},Int}()
      owner_p4est_gedge_subedge_to_hanging_edge = Dict{Tuple{Int,Int},Int}()
      num_hanging_faces[2] = 0
      ledge_to_cvertices = Gridap.ReferenceFEs.get_faces(HEX, 1, 0)
      # The following loop needs (1) the pairs to be traversed in increased order by cell ID;
      # (2) gridap cell vertices to be already completed
      for key in sort(collect(keys(hanging_edges_cell_ledge_to_owner_face_half)))
        (cell, ledge) = key
        (owner_p4est_gface_or_gedge, half) = hanging_edges_cell_ledge_to_owner_face_half[key]
        @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))] own_length=$(own_length(indices)) cell=$(cell) ledge=$(ledge) owner_p4est_gface_or_gedge=$(owner_p4est_gface_or_gedge) half=$(half)"
        if (half<0) # hanging edge is within a coarser face 
          owner_p4est_gface = owner_p4est_gface_or_gedge
          if !(haskey(owner_p4est_gface_half_to_hanging_edge, (owner_p4est_gface,half)))
            num_hanging_faces[2] += 1
            owner_p4est_gface_half_to_hanging_edge[(owner_p4est_gface,half)] = num_hanging_faces[2]
            if (!is_ghost(cell) || (haskey(regular_faces_p4est_to_gridap,owner_p4est_gface)))
              (owner_cell, p4est_lface) = p4est_gface_to_gcell_p4est_lface[owner_p4est_gface]
              push!(hanging_edges_owner_cell_and_lface,
                    (owner_cell, n_cell_vertices+n_cell_edges+PXEST_2_GRIDAP_FACE[p4est_lface],half))
            else
              push!(hanging_edges_owner_cell_and_lface,(-1, -1, -1))
            end 
          end
          start_gridap_edges = (cell - 1) * n_cell_edges
          gridap_cell_edges_data[start_gridap_edges+ledge] = num_regular_faces[2] + 
                        owner_p4est_gface_half_to_hanging_edge[(owner_p4est_gface,half)]

        else # hanging edge is within a coarser edge
          @assert half==1 || half==2
          owner_p4est_gedge = owner_p4est_gface_or_gedge
          owner_gedge_pair = (owner_p4est_gedge,half)
          if (haskey(owner_edge_subedge_to_cell_ledge,owner_gedge_pair))
            (owner_cell, owner_cell_ledge) = owner_edge_subedge_to_cell_ledge[owner_gedge_pair]
            if (owner_cell==cell)
              @assert owner_cell_ledge == ledge
              num_hanging_faces[2] += 1
              owner_p4est_gedge_subedge_to_hanging_edge[(owner_p4est_gedge,half)] = num_hanging_faces[2]
              if (!is_ghost(cell) || (haskey(regular_edges_p4est_to_gridap,owner_p4est_gedge)))
                (owner_cell, p4est_ledge) = p4est_gedge_to_gcell_p4est_ledge[owner_p4est_gedge]
                push!(hanging_edges_owner_cell_and_lface,
                     (owner_cell, n_cell_vertices+PXEST_2_GRIDAP_EDGE[p4est_ledge],half))
              else
                push!(hanging_edges_owner_cell_and_lface,(-1, -1, -1))
              end
              start_gridap_edges = (cell - 1) * n_cell_edges
              gridap_cell_edges_data[start_gridap_edges+ledge] = num_regular_faces[2] + num_hanging_faces[2] 
            else
              haskey_first_subedge  = haskey(owner_p4est_gedge_subedge_to_hanging_edge,(owner_p4est_gedge,1))
              haskey_second_subedge = haskey(owner_p4est_gedge_subedge_to_hanging_edge,(owner_p4est_gedge,2))
              if (!(is_ghost(cell)))
                @assert haskey_first_subedge && haskey_second_subedge
              else
                @assert haskey_first_subedge || haskey_second_subedge
              end 
              if (haskey_first_subedge && haskey_second_subedge)
                # The following code is required as we may have edges 
                # with different orientations at the inter-octree boundaries
                match=true
                start_gridap_vertices_cell = (cell - 1) * n_cell_vertices
                start_gridap_vertices_cell_owner = (owner_cell - 1) * n_cell_vertices
                for lvertex_cell in ledge_to_cvertices[ledge]
                  vertex_cell=gridap_cell_vertices_data[start_gridap_vertices_cell+lvertex_cell]
                  found=false
                  # Go over vertices of owner_cell_ledge in owner_cell 
                  for lvertex_owner_cell in ledge_to_cvertices[owner_cell_ledge]
                    vertex_owner_cell=gridap_cell_vertices_data[start_gridap_vertices_cell_owner+lvertex_owner_cell]
                    if (vertex_owner_cell==vertex_cell)
                      found=true
                      break
                    end
                  end
                  if (!found)
                    match=false
                    break
                  end
                end
                if (match)
                  owner_half=half
                else 
                  owner_half=half==1 ? 2 : 1
                end 
              elseif (haskey_first_subedge)
                owner_half=1
              elseif (haskey_second_subedge)
                owner_half=2
              end 
              @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))] cell=$(cell) ledge=$(ledge) owner_p4est_gface_or_gedge=$(owner_p4est_gface_or_gedge) half=$(half) owner_half=$(owner_half)"
              start_gridap_edges = (cell - 1) * n_cell_edges
              gridap_cell_edges_data[start_gridap_edges+ledge] = 
                 num_regular_faces[2] + owner_p4est_gedge_subedge_to_hanging_edge[(owner_p4est_gedge,owner_half)] 
            end
          end
        end
      end
      @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))] gridap_cell_edges_data: $(gridap_cell_edges_data)"
    end


    gridap_cell_faces = Vector{JaggedArray}(undef,Dc)
    gridap_cell_faces[1] = JaggedArray(gridap_cell_vertices_data,gridap_cell_vertices_ptrs)
    if (Dc==3)
      gridap_cell_faces[2] = JaggedArray(gridap_cell_edges_data,gridap_cell_edges_ptrs)
    end   
    gridap_cell_faces[Dc] = JaggedArray(gridap_cell_faces_data,gridap_cell_faces_ptrs)

    hanging_faces_glue      = Vector{Vector{Tuple}}(undef,Dc)
    hanging_faces_glue[1]   = hanging_vertices_owner_cell_and_lface
    if (Dc==3)
      hanging_faces_glue[2] = hanging_edges_owner_cell_and_lface
    end 
    hanging_faces_glue[Dc]  = hanging_faces_owner_cell_and_lface


    return num_regular_faces, 
           num_hanging_faces,
           gridap_cell_faces,
           hanging_faces_glue

  end |> tuple_of_arrays

  
  gridap_cell_faces_out  = Vector{MPIArray}(undef,Dc)
  for i=1:Dc
    gridap_cell_faces_out[i] = map(gridap_cell_faces) do  gridap_cell_faces
      gridap_cell_faces[i]
    end
  end
  non_conforming_glue=map(num_regular_faces,num_hanging_faces,hanging_faces_glue) do nrf, nhf, hfg
    NonConformingGlue(nrf, nhf, hfg)
  end 
  gridap_cell_faces_out,non_conforming_glue
 end


