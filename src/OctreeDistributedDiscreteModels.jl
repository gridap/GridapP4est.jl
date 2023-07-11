
const nothing_flag = Cint(0)
const refine_flag  = Cint(1)
const coarse_flag  = Cint(2)

mutable struct OctreeDistributedDiscreteModel{Dc,Dp,A,B,C,D,E} <: GridapDistributed.DistributedDiscreteModel{Dc,Dp}
  parts                       :: A
  dmodel                      :: B
  coarse_model                :: C
  ptr_pXest_connectivity      :: D
  ptr_pXest                   :: E

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
    C = typeof(coarse_model)
    D = typeof(ptr_pXest_connectivity)
    E = typeof(ptr_pXest)
    model = new{Dc,Dp,A,B,C,D,E}(parts,
                                 dmodel,
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
  coarse_model,
  ptr_pXest_connectivity,
  ptr_pXest,
  owns_ptr_pXest_connectivity,
  gc_ref) where {Dc,Dp}

  return OctreeDistributedDiscreteModel(Dc,
                                        Dp,
                                        parts,
                                        dmodel,
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
  if i_am_in(comm)
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

    return OctreeDistributedDiscreteModel(Dc,
                                          Dp,
                                          parts,
                                          dmodel,
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
    p4est_partition(ptr_pXest, 0, C_NULL)
  else
    @assert false
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

function pXest_refine!(::Type{Val{Dc}}, ptr_pXest, refine_fn_c) where Dc
  if (Dc==2)
    p4est_refine(ptr_pXest, Cint(0), refine_fn_c, C_NULL)
  else
    p8est_refine(ptr_pXest, Cint(0), refine_fn_c, C_NULL)
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
    @assert false
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


function _compute_fine_to_coarse_model_glue(
         cparts,
         cmodel::Union{Nothing,GridapDistributed.DistributedDiscreteModel{Dc}},
         fmodel::GridapDistributed.DistributedDiscreteModel{Dc}) where Dc

  # Fill data for owned (from coarse cells owned by the processor)
  fgids = get_cell_gids(fmodel)
  f1,f2,f3, cgids_snd, cgids_rcv = map(fmodel.models,
                                       partition(fgids)) do fmodel, fpartition
    if (!(i_am_in(cparts)))
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
  # dfcell_to_child_id = map(f3) do fcell_to_child_id
  #  !isa(fcell_to_child_id,Nothing) ? fcell_to_child_id : Int[]
  # end
  # cache=fetch_vector_ghost_values_cache(dfcell_to_child_id,partition(fgids))
  # fetch_vector_ghost_values!(dfcell_to_child_id,cache) |> wait
  
  # Note: Reversing snd and rcv 
  parts_rcv, parts_snd = assembly_neighbors(partition(fgids))
  PArrays.exchange_fetch!(cgids_rcv,cgids_snd,ExchangeGraph(parts_snd,parts_rcv))

  map(f1,cgids_rcv) do fine_to_coarse_faces_map, cgids_rcv
    if (i_am_in(cparts))
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
    if (!(i_am_in(cparts)))
      nothing
    else
      cmodel_local = PArrays.getany(cmodel.models)
      num_cells_coarse = num_cells(cmodel_local)
      reffe  = LagrangianRefFE(Float64,QUAD,1)
      rrules = Fill(Gridap.Adaptivity.RefinementRule(reffe,2),num_cells_coarse)
      #println("fcell_to_child_id: ", fcell_to_child_id)
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

function Gridap.Adaptivity.refine(model::OctreeDistributedDiscreteModel{Dc,Dp}, parts=nothing) where {Dc,Dp}
   old_comm = model.parts.comm
   if (i_am_in(old_comm))
     # Copy and refine input p4est
     ptr_new_pXest = pXest_copy(Val{Dc}, model.ptr_pXest)
     pXest_uniformly_refine!(Val{Dc}, ptr_new_pXest)
   else
     ptr_new_pXest = nothing
   end

   new_comm = isa(parts,Nothing) ? old_comm : parts.comm
   if i_am_in(new_comm)
      if !isa(parts,Nothing)
        aux = ptr_new_pXest
        ptr_new_pXest = _p4est_to_new_comm(ptr_new_pXest,
                                           model.ptr_pXest_connectivity,
                                           model.parts.comm,
                                           parts.comm)
        if i_am_in(old_comm)
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

      ref_model = OctreeDistributedDiscreteModel(Dc,Dp,
                                     new_parts,
                                     fmodel,
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

function Gridap.Adaptivity.refine(model::OctreeDistributedDiscreteModel{Dc,Dp}, 
                                  refinement_and_coarsening_flags::MPIArray{<:Vector}) where {Dc,Dp}

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
    end                                     

    map(model.dmodel.models,refinement_and_coarsening_flags) do lmodel, flags
      # The length of the local flags array has to match the number of 
      # cells in the model. This includes both owned and ghost cells. 
      # Only the flags for owned cells is actually taken into account. 
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
    
    # Copy and refine input p4est
    ptr_new_pXest = pXest_copy(Val{Dc}, model.ptr_pXest)
    pXest_refine!(Val{Dc}, ptr_new_pXest, refine_callback_2d_c)

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


     # To-DO: does not work OK in the non-conforming case 
     #  dglue = _compute_fine_to_coarse_model_glue(model.parts,
     #                                             model.dmodel,
     #                                             fmodel)
     ref_model = OctreeDistributedDiscreteModel(Dc,Dp,
                                    model.parts,
                                    fmodel,
                                    model.coarse_model,
                                    model.ptr_pXest_connectivity,
                                    ptr_new_pXest,
                                    false,
                                    model)
     return ref_model, non_conforming_glue
  # else
  #  new_parts = isa(parts,Nothing) ? model.parts : parts
  #  return VoidOctreeDistributedDiscreteModel(model,new_parts), nothing
  # end
end

function Gridap.Adaptivity.coarsen(model::OctreeDistributedDiscreteModel{Dc,Dp}) where {Dc,Dp}
  comm = model.parts.comm
  if (i_am_in(comm))
    # Copy and refine input p4est
    ptr_new_pXest = pXest_copy(Val{Dc}, model.ptr_pXest)
    pXest_uniformly_coarsen!(Val{Dc}, ptr_new_pXest)
  else
    ptr_new_pXest=nothing
  end

  if (i_am_in(comm))
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

     c_octree_model = OctreeDistributedDiscreteModel(Dc,Dp,
                                    model.parts,
                                    cmodel,
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
  if (i_am_in(new_comm))
    new_comm_num_parts    = num_parts(new_comm)
    global_first_quadrant = Vector{P4est_wrapper.p4est_gloidx_t}(undef,new_comm_num_parts+1)

    pXest_conn = ptr_pXest_conn[]
    pertree = Vector{P4est_wrapper.p4est_gloidx_t}(undef,pXest_conn.num_trees+1)
    if (i_am_in(old_comm))
      pXest = ptr_pXest[]
      old_comm_num_parts = num_parts(old_comm)
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
  @assert i_am_in(old_comm)
  pXest = ptr_pXest[]
  pXest_conn = ptr_pXest_conn[]

  pertree = Vector{P4est_wrapper.p4est_gloidx_t}(undef,pXest_conn.num_trees+1)
  p4est_comm_count_pertree(ptr_pXest,pertree)

  if (i_am_in(new_comm))
    new_comm_num_parts = num_parts(new_comm)
    global_first_quadrant = Vector{P4est_wrapper.p4est_gloidx_t}(undef,new_comm_num_parts+1)
    pXest=ptr_pXest[]
    old_comm_num_parts = num_parts(old_comm)
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


function _p4est_tree_array_index(::Type{Val{Dc}},trees,itree) where Dc
  if (Dc==2)
    return p4est_tree_array_index(trees,itree)
  elseif (Dc==3)
    return p8est_tree_array_index(trees,itree)
  end
end

function _p4est_comm_find_owner(::Type{Val{Dc}},ptr_pXest,itree,quad,guess) where Dc
  if (Dc==2)
    return p4est_comm_find_owner(ptr_pXest,itree,quad,guess)
  elseif (Dc==3)
    return p8est_comm_find_owner(ptr_pXest,itree,quad,guess)
  end
end

function _p4est_quadrant_array_index(::Type{Val{Dc}}, quadrants, iquad) where Dc
  if (Dc==2)
    return p4est_quadrant_array_index(quadrants, iquad)
  elseif (Dc==3)
    return p8est_quadrant_array_index(quadrants, iquad)
  end
end


function _p4est_quadrant_is_equal(::Type{Val{Dc}},  q1, q2) where Dc
  if (Dc==2)
    return p4est_quadrant_is_equal(q1,q2)
  elseif (Dc==3)
    return p8est_quadrant_is_equal(q1, q2)
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
    tree = _p4est_tree_array_index(Val{Dc},pXest_old.trees,itree)[]
    num_quads = Cint(tree.quadrants.elem_count)

    for iquad = 0:num_quads-1
      q = _p4est_quadrant_array_index(Val{Dc},tree.quadrants, iquad)
      new_rank = _p4est_comm_find_owner(Val{Dc},ptr_pXest_new,itree,q,0)

      if (new_rank != my_rank)
        if (!(new_rank+1 in keys(ranks_count)))
          push!(lst_ranks,new_rank+1)
          ranks_count[new_rank+1] = 0
        end
        ranks_count[new_rank+1] += 1
        old2new[current_old_quad_index] = 0
      else
        current_new_quad_index = 1
        new_tree = _p4est_tree_array_index(Val{Dc},pXest_new.trees,pXest_new.first_local_tree)[]
        for t = pXest_new.first_local_tree:pXest_new.last_local_tree
          new_tree = _p4est_tree_array_index(Val{Dc},pXest_new.trees,t)[]
          if t == itree
            break
          end
          current_new_quad_index += Cint(new_tree.quadrants.elem_count)
        end
        found = false
        num_quads_new = Cint(new_tree.quadrants.elem_count)
        for iquad_new = 0:num_quads_new-1
          q_new = _p4est_quadrant_array_index(Val{Dc},new_tree.quadrants, iquad_new)
          found = _p4est_quadrant_is_equal(Val{Dc},q,q_new)!=0
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
  @assert i_am_in(partsA.comm) || i_am_in(partsB.comm)
  is_included(partsA.comm,partsB.comm)
end

function is_included(commA::MPI.Comm,commB::MPI.Comm)
  @assert i_am_in(commA) || i_am_in(commB)
  num_partsA=num_parts(commA)
  num_partsB=num_parts(commB)
  if (num_partsA==num_partsB)
    return false
  end
  if (i_am_in(commA) && i_am_in(commB))
    result=num_partsA < num_partsB
    if (result)
      result=MPI.Allreduce(Int8(result),MPI.LOR,commB)
    else
      result=MPI.Allreduce(Int8(result),MPI.LOR,commA)
    end
  else
    result=false
    if (i_am_in(commB))
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
function GridapDistributed.redistribute(model::OctreeDistributedDiscreteModel{Dc,Dp}, parts_redistributed_model=model.parts) where {Dc,Dp}
  parts = (parts_redistributed_model === model.parts) ? model.parts : parts_redistributed_model
  comm  = parts.comm
  if (i_am_in(model.parts.comm) || i_am_in(parts.comm))
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

  glue = GridapDistributed.RedistributeGlue(parts_rcv,parts_snd,lids_rcv,lids_snd,old2new,new2old)

  # Extract ghost and lnodes
  ptr_pXest_ghost  = setup_pXest_ghost(Val{Dc}, ptr_pXest_new)
  ptr_pXest_lnodes = setup_pXest_lnodes(Val{Dc}, ptr_pXest_new, ptr_pXest_ghost)

  # Build fine-grid mesh
  fmodel = setup_distributed_discrete_model(Val{Dc},
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
  if (i_am_in(subset_comm))
    # This piece of code replicates the logic behind the
    # "p4est_partition_cut_gloidx" function in the p4est library
    psub=PArrays.getany(parts_redistributed_model)
    Psub=num_parts(subset_comm)
    first_global_quadrant=Int64((Float64(N)*Float64(psub-1))/(Float64(Psub)))
    @assert first_global_quadrant>=0 && first_global_quadrant<N
    # print("$(N) $(Psub) $(psub): $(first_global_quadrant)","\n")
  else
    first_global_quadrant=N
  end

  Psup=num_parts(supset_comm)
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

  glue = GridapDistributed.RedistributeGlue(parts_rcv,parts_snd,lids_rcv,lids_snd,old2new,new2old)

  if (i_am_in(subset_comm))
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

    red_model = OctreeDistributedDiscreteModel(Dc,Dp,
                                               parts_redistributed_model,
                                               fmodel,
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
        hanging_edge[e+1] = (hanging_edge[e] == -1) ? 0 : 2
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

  num_regular_faces,
  num_hanging_faces,
  gridap_cell_faces,
  hanging_faces_glue =
    generate_cell_faces_and_non_conforming_glue(Val{Dc},ptr_pXest_lnodes, cell_prange)

  println("### faces ###")
  println("num_regular_faces: $(num_regular_faces)")
  println("num_hanging_faces: $(num_hanging_faces)")
  println("gridap_cell_faces: $(gridap_cell_faces)")
  println(hanging_faces_glue)

  nlvertices = map(num_regular_faces[1],num_hanging_faces[1]) do nrv,nhv
    nrv+nhv
  end

  node_coordinates=generate_node_coordinates(Val{Dc},
                                             gridap_cell_faces[1],
                                             nlvertices,
                                             ptr_pXest_connectivity,
                                             ptr_pXest,
                                             ptr_pXest_ghost)


  println("node_coordinates: $(node_coordinates)")

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

  _set_hanging_labels!(face_labeling,num_regular_faces,num_hanging_faces)

  discretemodel=map(grid,topology,face_labeling) do grid, topology, face_labeling
    Gridap.Geometry.UnstructuredDiscreteModel(grid,topology,face_labeling)
  end
  GridapDistributed.DistributedDiscreteModel(discretemodel,cell_prange), 
    (num_regular_faces,num_hanging_faces,gridap_cell_faces,hanging_faces_glue)
end

function _set_hanging_labels!(face_labeling,num_regular_faces,num_hanging_faces)
  max_entity_ids = map(face_labeling) do face_labeling
    max_entity_id = typemin(eltype(first(face_labeling.d_to_dface_to_entity))) 
    for i=1:length(face_labeling.d_to_dface_to_entity)
      max_entity_id=max(maximum(face_labeling.d_to_dface_to_entity[i]),max_entity_id)
    end
    max_entity_id
  end
  max_entity_id = MPI.Allreduce(max_entity_ids.part,MPI.MAX,max_entity_ids.comm)
  
  hanging_entitity_ids = Dict{Int,Bool}()
  for i=1:length(num_hanging_faces)
     map(face_labeling,
               num_regular_faces[i],
               num_hanging_faces[i]) do face_labeling, num_regular_faces, num_hanging_faces
      for j=num_regular_faces+1:num_regular_faces+num_hanging_faces
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

function _build_map_from_faces_to_cell_lface(lnodes)
  n_cell_faces    = num_cell_faces(Val{2})

  element_nodes = unsafe_wrap(Array, lnodes.element_nodes, lnodes.vnodes * lnodes.num_local_elements)
  face_code = unsafe_wrap(Array, lnodes.face_code, lnodes.num_local_elements)
  hanging_face = Vector{Cint}(undef, n_cell_faces)

  # Build a map from faces to (cell,lface)
  p4est_gface_to_gcell_p4est_lface = Dict{Int,Tuple{Int,Int}}()
  for cell = 1:lnodes.num_local_elements
    start = (cell - 1) * lnodes.vnodes + 1
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

function _build_map_from_faces_edges_to_cell_lface_ledge(lnodes)
    n_cell_faces    = num_cell_faces(Val{3})
    n_cell_edges    = num_cell_edges(Val{3})

    element_nodes = unsafe_wrap(Array, lnodes.element_nodes, lnodes.vnodes * lnodes.num_local_elements)
    face_code = unsafe_wrap(Array, lnodes.face_code, lnodes.num_local_elements)
    hanging_face = Vector{Cint}(undef, n_cell_faces)
    hanging_edge = Vector{Cint}(undef, n_cell_edges)

    # Build a map from faces to (cell,lface)
    p4est_gface_to_gcell_p4est_lface = Dict{Int,Tuple{Int,Int}}()
    p4est_gedge_to_gcell_p4est_ledge = Dict{Int,Tuple{Int,Int}}()
    for cell = 1:lnodes.num_local_elements
      start = (cell - 1) * lnodes.vnodes + 1
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
  hanging_faces_glue = map(cell_prange.partition) do indices

    num_regular_faces = Vector{Int}(undef, Dc)
    num_hanging_faces = Vector{Int}(undef, Dc)

    # Count regular vertices
    num_regular_faces[1] = 0
    regular_vertices_p4est_to_gridap = Dict{Int,Int}()

    num_regular_faces[Dc] = 0
    regular_faces_p4est_to_gridap = Dict{Int,Int}()

    if (Dc==2)
      p4est_gface_to_gcell_p4est_lface = _build_map_from_faces_to_cell_lface(lnodes)
    else 
      p4est_gface_to_gcell_p4est_lface, 
         p4est_gedge_to_gcell_p4est_ledge = _build_map_from_faces_edges_to_cell_lface_ledge(lnodes)
    end 

    PXEST_2_GRIDAP_VERTEX = pXest_2_gridap_vertex(Val{Dc})
    PXEST_2_GRIDAP_FACE   = pXest_2_gridap_facet(Val{Dc})
    PXEST_2_GRIDAP_EDGE   = p8est_2_gridap_edge()

    n = length(indices.lid_to_part)
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

    gridap_cell_vertices_data = Vector{Int}(undef, lnodes.num_local_elements * n_cell_vertices)
    gridap_cell_vertices_data .= -1

    gridap_cell_faces_data = Vector{Int}(undef, lnodes.num_local_elements * n_cell_faces)
    gridap_cell_faces_data .= -1

    if (Dc==3)
      num_regular_faces[2] = 0

      gridap_cell_edges_ptrs = Vector{Int32}(undef,n+1)
      gridap_cell_edges_ptrs[1]=1
      for i=1:n
        gridap_cell_edges_ptrs[i+1]=gridap_cell_edges_ptrs[i]+n_cell_edges
      end
      gridap_cell_edges_data = Vector{Int}(undef, lnodes.num_local_elements * n_cell_edges)
      gridap_cell_edges_data .= -1
      hanging_edges_pairs_to_owner_face = Dict{Tuple{Int,Int},Tuple{Int,Int}}()
      hanging_vertices_pairs_to_owner_edge = Dict{Tuple{Int,Int},Int}()
      regular_edges_p4est_to_gridap = Dict{Int,Int}()
    end 

    for cell = 1:lnodes.num_local_elements
      start                 = (cell - 1) * lnodes.vnodes + 1
      start_gridap_vertices = (cell - 1) * n_cell_vertices
      start_gridap_faces    = (cell - 1) * n_cell_faces

      p4est_cell_faces = view(element_nodes, start:start+n_cell_faces-1)
      p4est_cell_vertices = view(element_nodes, 
                                 start+n_cell_faces+n_cell_edges:start+n_cell_faces+n_cell_edges+n_cell_vertices-1)

      gridap_cell_vertices = view(gridap_cell_vertices_data,
        start_gridap_vertices+1:start_gridap_vertices+n_cell_vertices)
      gridap_cell_faces = view(gridap_cell_faces_data,
        start_gridap_faces+1:start_gridap_faces+n_cell_faces)

      if (Dc==2)  
        has_hanging = p4est_lnodes_decode(face_code[cell], hanging_face)
      else
        has_hanging = p8est_lnodes_decode(face_code[cell], hanging_face, hanging_edge)
        start_gridap_edges = (cell-1)*n_cell_edges
        gridap_cell_edges = view(gridap_cell_edges_data, start_gridap_edges+1:start_gridap_edges+n_cell_edges)
        p4est_cell_edges = view(element_nodes, 
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
                # Identify the two edges which are hanging within the face
                hanging_edges_pairs_to_owner_face[(cell, PXEST_2_GRIDAP_EDGE[p4est_ledge+1])] =
                    (owner_face,-p8est_subface_to_hanging_edges_within_face[half+1,i])
                gridap_cell_edges[PXEST_2_GRIDAP_EDGE[p4est_ledge+1]] = hanging_edge_from_face_code
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
                owner_edge = p4est_cell_edges[p4est_ledge]
                hanging_vertices_pairs_to_owner_edge[(cell, 
                                                      PXEST_2_GRIDAP_VERTEX[p4est_hanging_lvertex+1])] = owner_edge

                # Process hanging edge
                hanging_edges_pairs_to_owner_face[(cell, PXEST_2_GRIDAP_EDGE[p4est_ledge])] = 
                   (owner_edge,regular_vertex_lvertex_within_edge+1)
              end 
            end 
        end
      end
    end
  end 

    # Go over all touched hanging faces and start 
    # assigning IDs from the last num_regular_faces ID
    # For each hanging face, keep track of (owner_cell,lface)
    hanging_faces_owner_cell_and_lface =
      Vector{Tuple{Int,Int,Int}}(undef, length(keys(hanging_faces_pairs_to_owner_face)))
    num_hanging_faces[Dc] = 0
    for key in keys(hanging_faces_pairs_to_owner_face)
      (cell, lface) = key
      (owner_p4est_gface, half) = hanging_faces_pairs_to_owner_face[key]
      owner_gridap_gface = regular_faces_p4est_to_gridap[owner_p4est_gface]
      num_hanging_faces[Dc] += 1
      start_gridap_faces = (cell - 1) * n_cell_faces
      gridap_cell_faces_data[start_gridap_faces+lface] = num_regular_faces[Dc] + num_hanging_faces[Dc]
      (owner_cell, p4est_lface) = p4est_gface_to_gcell_p4est_lface[owner_p4est_gface]
      hanging_faces_owner_cell_and_lface[num_hanging_faces[Dc]] =
        (owner_cell, n_cell_vertices+n_cell_edges+PXEST_2_GRIDAP_FACE[p4est_lface], half)
    end

    println("AAA: $(gridap_cell_faces_data)")

    # Go over all touched hanging vertices and start 
    # assigning IDs from the last num_regular_vertices ID
    # For each hanging vertex, keep track of (owner_cell,lface)
    num_hanging_faces[1] = 0
    hanging_vertices_owner_cell_and_lface = Tuple{Int,Int,Int}[]
    owner_gridap_gface_to_hanging_vertex = Dict{Int,Int}()
    half=1
    for key in keys(hanging_vertices_pairs_to_owner_face)
      (cell, lvertex) = key
      owner_p4est_gface = hanging_vertices_pairs_to_owner_face[key]
      owner_gridap_gface = regular_faces_p4est_to_gridap[owner_p4est_gface]
      if !(haskey(owner_gridap_gface_to_hanging_vertex, owner_gridap_gface))
        num_hanging_faces[1] += 1
        owner_gridap_gface_to_hanging_vertex[owner_gridap_gface] = num_hanging_faces[1]
        (owner_cell, p4est_lface) = p4est_gface_to_gcell_p4est_lface[owner_p4est_gface]
        push!(hanging_vertices_owner_cell_and_lface,
          (owner_cell, n_cell_vertices+n_cell_edges+PXEST_2_GRIDAP_FACE[p4est_lface],half))
      end
      start_gridap_vertices = (cell - 1) * n_cell_vertices
      gridap_cell_vertices_data[start_gridap_vertices+lvertex] = num_regular_faces[1] +
                                                             owner_gridap_gface_to_hanging_vertex[owner_gridap_gface]
    end

    if (Dc==3)
      # Go over all touched hanging edges and start 
      # assigning IDs from the last num_regular_edge ID
      # For each hanging edge, keep track of (owner_cell,lface/ledge)
      hanging_edges_owner_cell_and_lface = Tuple{Int,Int,Int}[]

      owner_gridap_gface_half_to_hanging_edge = Dict{Tuple{Int,Int},Int}()
      owner_gridap_gedge_half_to_hanging_edge = Dict{Tuple{Int,Int},Int}()
      num_hanging_faces[2] = 0
      for key in keys(hanging_edges_pairs_to_owner_face)
        (cell, ledge) = key
        (owner_p4est_gface_or_gedge, half) = hanging_edges_pairs_to_owner_face[key]
        if (half<0) # hanging edge is within a coarser face 
          owner_gridap_gface = regular_faces_p4est_to_gridap[owner_p4est_gface_or_gedge]
          if !(haskey(owner_gridap_gface_half_to_hanging_edge, (owner_gridap_gface,half)))
            num_hanging_faces[2] += 1
            owner_gridap_gface_half_to_hanging_edge[(owner_gridap_gface,half)] = num_hanging_faces[2]
            (owner_cell, p4est_lface) = p4est_gface_to_gcell_p4est_lface[owner_p4est_gface_or_gedge]
            push!(hanging_edges_owner_cell_and_lface,
              (owner_cell, n_cell_vertices+n_cell_edges+PXEST_2_GRIDAP_FACE[p4est_lface],half))
          end
          start_gridap_edges = (cell - 1) * n_cell_edges
          gridap_cell_edges_data[start_gridap_edges+ledge] = num_regular_faces[2] + 
                        owner_gridap_gface_half_to_hanging_edge[(owner_gridap_gface,half)]

        else          # hanging edge is within a coarser edge
          @assert half==1 || half==2
          owner_gridap_gedge = regular_edges_p4est_to_gridap[owner_p4est_gface_or_gedge]
          if !(haskey(owner_gridap_gedge_half_to_hanging_edge, (owner_gridap_gedge,half)))
            num_hanging_faces[2] += 1
            owner_gridap_gedge_half_to_hanging_edge[(owner_gridap_gedge,half)] = num_hanging_faces[2]
            (owner_cell, p4est_ledge) = p4est_gedge_to_gcell_p4est_ledge[owner_p4est_gface_or_gedge]
            push!(hanging_edges_owner_cell_and_lface,
              (owner_cell, n_cell_vertices+PXEST_2_GRIDAP_EDGE[p4est_ledge],half))
          end
          start_gridap_edges = (cell - 1) * n_cell_edges
          gridap_cell_edges_data[start_gridap_edges+ledge] = num_regular_faces[2] + 
               owner_gridap_gedge_half_to_hanging_edge[(owner_gridap_gedge,half)]
        end
      end 
      
      half=1
      owner_gridap_gedge_to_hanging_vertex = Dict{Int,Int}()
      for key in keys(hanging_vertices_pairs_to_owner_edge)
        (cell, lvertex) = key
        owner_p4est_gedge = hanging_vertices_pairs_to_owner_edge[key]
        owner_gridap_gedge = regular_edges_p4est_to_gridap[owner_p4est_gedge]
        if !(haskey(owner_gridap_gedge_to_hanging_vertex, owner_gridap_gedge))
          num_hanging_faces[1] += 1
          owner_gridap_gedge_to_hanging_vertex[owner_gridap_gedge] = num_hanging_faces[1]
          (owner_cell, p4est_ledge) = p4est_gedge_to_gcell_p4est_ledge[owner_p4est_gedge]
          push!(hanging_vertices_owner_cell_and_lface,
            (owner_cell, n_cell_vertices+PXEST_2_GRIDAP_EDGE[p4est_ledge],half))
        end
        start_gridap_vertices = (cell - 1) * n_cell_vertices
        gridap_cell_vertices_data[start_gridap_vertices+lvertex] = num_regular_faces[1] +
                                                       owner_gridap_gedge_to_hanging_vertex[owner_gridap_gedge]
      end
    end

    println("XXX: $(gridap_cell_vertices_data)")

    gridap_cell_faces = Vector{Table}(undef,Dc)
    gridap_cell_faces[1] = Table(gridap_cell_vertices_data,gridap_cell_vertices_ptrs)
    if (Dc==3)
      gridap_cell_faces[2] = Table(gridap_cell_edges_data,gridap_cell_edges_ptrs)
    end   
    gridap_cell_faces[Dc] = Table(gridap_cell_faces_data,gridap_cell_faces_ptrs)

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

  end

  num_regular_faces_out  = Vector{MPIArray}(undef,Dc)
  num_hanging_faces_out  = Vector{MPIArray}(undef,Dc)
  gridap_cell_faces_out  = Vector{MPIArray}(undef,Dc)
  hanging_faces_glue_out = Vector{MPIArray}(undef,Dc)

  for i=1:Dc 
    num_regular_faces_out[i] = map(num_regular_faces) do  num_regular_faces
      num_regular_faces[i]
    end
    num_hanging_faces_out[i] = map(num_hanging_faces) do  num_hanging_faces
      num_hanging_faces[i]
    end
    gridap_cell_faces_out[i] = map(gridap_cell_faces) do  gridap_cell_faces
      gridap_cell_faces[i]
    end
    hanging_faces_glue_out[i] = map(hanging_faces_glue) do hanging_faces_glue
      hanging_faces_glue[i]
    end
  end
  num_regular_faces_out, 
  num_hanging_faces_out, 
  gridap_cell_faces_out, 
  hanging_faces_glue_out
 end


