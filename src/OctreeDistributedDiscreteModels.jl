

mutable struct OctreeDistributedDiscreteModel{Dc,Dp,A,B,C,D,E} <: GridapDistributed.AbstractDistributedDiscreteModel{Dc,Dp}
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
    dmodel::Union{GridapDistributed.AbstractDistributedDiscreteModel,Nothing},
    coarse_model,
    ptr_pXest_connectivity,
    ptr_pXest,
    owns_ptr_pXest_connectivity::Bool,
    gc_ref)

    if (isa(dmodel,GridapDistributed.AbstractDistributedDiscreteModel))
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
  dmodel::GridapDistributed.AbstractDistributedDiscreteModel{Dc,Dp},
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


function OctreeDistributedDiscreteModel(parts::MPIData{<:Integer},
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
    parts::MPIData{<:Integer},
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

# AbstractDistributedDiscreteModel API implementation

Gridap.Geometry.num_cells(model::OctreeDistributedDiscreteModel) = Gridap.Geometry.num_cells(model.dmodel)
Gridap.Geometry.num_facets(model::OctreeDistributedDiscreteModel) = Gridap.Geometry.num_facets(model.dmodel)
Gridap.Geometry.num_edges(model::OctreeDistributedDiscreteModel) = Gridap.Geometry.num_edges(model.dmodel)
Gridap.Geometry.num_vertices(model::OctreeDistributedDiscreteModel) = Gridap.Geometry.num_vertices(model.dmodel)
Gridap.Geometry.num_faces(model::OctreeDistributedDiscreteModel) = Gridap.Geometry.num_faces(model.dmodel)
Gridap.Geometry.get_grid(model::OctreeDistributedDiscreteModel) = Gridap.Geometry.get_grid(model.dmodel)
Gridap.Geometry.get_grid_topology(model::OctreeDistributedDiscreteModel) = Gridap.Geometry.get_grid_topology(model.dmodel)
Gridap.Geometry.get_face_labeling(model::OctreeDistributedDiscreteModel) = Gridap.Geometry.get_face_labeling(model.dmodel)

GridapDistributed.get_parts(model::OctreeDistributedDiscreteModel) = model.parts
GridapDistributed.local_views(model::OctreeDistributedDiscreteModel) = GridapDistributed.local_views(model.dmodel)
GridapDistributed.get_cell_gids(model::OctreeDistributedDiscreteModel) = GridapDistributed.get_cell_gids(model.dmodel)
GridapDistributed.get_face_gids(model::OctreeDistributedDiscreteModel,dim::Integer) = GridapDistributed.get_face_gids(model.dmodel,dim)
GridapDistributed.generate_gids(model::OctreeDistributedDiscreteModel,spaces) = GridapDistributed.generate_gids(model.dmodel,spaces)

# Garbage collection

function octree_distributed_discrete_model_free!(model::VoidOctreeDistributedDiscreteModel{Dc}) where Dc
  # parts = get_parts(model)
  # if i_am_in(parts)
    if (model.owns_ptr_pXest_connectivity)
      pXest_connectivity_destroy(Val{Dc},model.ptr_pXest_connectivity)
    end
  # end
  return nothing
end

function octree_distributed_discrete_model_free!(model::OctreeDistributedDiscreteModel{Dc}) where Dc
  parts = get_parts(model)
  if i_am_in(parts)
    pXest_destroy(Val{Dc},model.ptr_pXest)
    if (model.owns_ptr_pXest_connectivity)
      pXest_connectivity_destroy(Val{Dc},model.ptr_pXest_connectivity)
    end
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
    @assert false
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

function pXest_refine!(::Type{Val{Dc}}, ptr_pXest) where Dc
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

function pXest_coarsen!(::Type{Val{Dc}}, ptr_pXest) where Dc
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
  f1,f2,f3, cgids_snd, cgids_rcv = map_parts(fmodel.models,fgids.partition) do fmodel, fpartition
    if (!(i_am_in(cparts)))
      # cmodel might be distributed among less processes than fmodel
      nothing, nothing, nothing, Int[], Int[]
    else
      cgids        = get_cell_gids(cmodel)
      cmodel_local = cmodel.models.part
      cpartition   = cgids.partition.part
      fine_to_coarse_faces_map,
        fine_to_coarse_faces_dim,
          fcell_to_child_id =
          _process_owned_cells_fine_to_coarse_model_glue(cmodel_local,fmodel,cpartition,fpartition)

      lids_snd    = fgids.exchanger.lids_snd.part
      lids_rcv    = fgids.exchanger.lids_rcv.part
      cgids_data  = cpartition.lid_to_gid[fine_to_coarse_faces_map[Dc+1][lids_snd.data]]
      cgids_snd   = Table(cgids_data,lids_snd.ptrs)
      cgids_rcv   = Table(Vector{Int}(undef,length(lids_rcv.data)),lids_rcv.ptrs)

      fine_to_coarse_faces_map, fine_to_coarse_faces_dim, fcell_to_child_id, cgids_snd, cgids_rcv
    end
  end

  # Nearest Neighbors comm: Get data for ghosts (from coarse cells owned by neighboring processors)
  dfcell_to_child_id = map_parts(f3) do fcell_to_child_id
    !isa(fcell_to_child_id,Nothing) ? fcell_to_child_id : Int[]
  end
  exchange!(dfcell_to_child_id, fgids.exchanger)
  tout = PArrays.async_exchange!(cgids_rcv,
                         cgids_snd,
                         fgids.exchanger.parts_rcv,
                         fgids.exchanger.parts_snd)
  map_parts(schedule,tout)
  map_parts(wait,tout)

  map_parts(f1,cgids_rcv) do fine_to_coarse_faces_map, cgids_rcv
    if (i_am_in(cparts))
      cgids      = get_cell_gids(cmodel)
      cpartition = cgids.partition.part
      lids_rcv   = fgids.exchanger.lids_rcv.part

      for i in 1:length(lids_rcv.data)
        lid = lids_rcv.data[i]
        gid = cgids_rcv.data[i]
        fine_to_coarse_faces_map[Dc+1][lid] = cpartition.gid_to_lid[gid]
      end
    end
  end

  # Create distributed glue
  map_parts(f1,f2,f3) do fine_to_coarse_faces_map, fine_to_coarse_faces_dim, fcell_to_child_id
    if (!(i_am_in(cparts)))
      nothing
    else
      cmodel_local = cmodel.models.part
      num_cells_coarse = num_cells(cmodel_local)
      reffe  = LagrangianRefFE(Float64,QUAD,1)
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

  num_f_cells   = num_cells(fmodel)             # Number of fine cells (owned+ghost)
  num_o_c_cells = length(cpartition.oid_to_lid) # Number of coarse cells (owned)

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
     pXest_refine!(Val{Dc}, ptr_new_pXest)
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

function Gridap.Adaptivity.coarsen(model::OctreeDistributedDiscreteModel{Dc,Dp}) where {Dc,Dp}
  comm = model.parts.comm
  if (i_am_in(comm))
    # Copy and refine input p4est
    ptr_new_pXest = pXest_copy(Val{Dc}, model.ptr_pXest)
    pXest_coarsen!(Val{Dc}, ptr_new_pXest)
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

    #print("$(global_first_quadrant) " * " \n")
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

  lst_ranks,PartitionedArrays.Table(local_ids,ptr_ranks),old2new
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

function _redistribute_parts_supset_parts_redistributed(model::OctreeDistributedDiscreteModel{Dc,Dp}, parts_redistributed_model) where {Dc,Dp}
  @assert model.parts !== parts_redistributed_model

  subset_comm = parts_redistributed_model.comm
  supset_comm = model.parts.comm
  N=num_cells(model)
  if (i_am_in(subset_comm))
    # This piece of code replicates the logic behind the
    # "p4est_partition_cut_gloidx" function in the p4est library
    psub=parts_redistributed_model.part
    Psub=num_parts(subset_comm)
    first_global_quadrant=Int64((Float64(N)*Float64(psub-1))/(Float64(Psub)))
    @assert first_global_quadrant>=0 && first_global_quadrant<N
    # print("$(N) $(Psub) $(psub): $(first_global_quadrant)","\n")
  else
    first_global_quadrant=N
  end

  Psup=num_parts(supset_comm)
  psub=model.parts.part
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
  lids_rcv, parts_rcv = map_parts(parts) do _
    lids_rcv, parts_rcv
  end
  lids_snd, parts_snd = map_parts(parts) do _
    lids_snd, parts_snd
  end
  old2new, new2old = map_parts(parts) do _
    old2new,new2old
  end
  lids_rcv, parts_rcv, lids_snd, parts_snd, old2new, new2old
end
