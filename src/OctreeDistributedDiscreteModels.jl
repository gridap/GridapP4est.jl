

struct OctreeDistributedDiscreteModel{Dc,Dp,A,B,C,D,E} <: GridapDistributed.AbstractDistributedDiscreteModel{Dc,Dp}
  parts                  :: A
  dmodel                 :: B
  coarse_model           :: C
  ptr_pXest_connectivity :: D
  ptr_pXest              :: E

  function OctreeDistributedDiscreteModel(
    Dc::Int, 
    Dp::Int,
    parts,
    dmodel::Union{GridapDistributed.AbstractDistributedDiscreteModel,Nothing},
    coarse_model,
    ptr_pXest_connectivity,
    ptr_pXest)

    if (isa(dmodel,GridapDistributed.AbstractDistributedDiscreteModel))
      Gridap.Helpers.@check Dc == Gridap.Geometry.num_cell_dims(dmodel)
      Gridap.Helpers.@check Dc == Gridap.Geometry.num_point_dims(dmodel)
    end
  
    A = typeof(parts)
    B = typeof(dmodel)
    C = typeof(coarse_model)
    D = typeof(ptr_pXest_connectivity)
    E = typeof(ptr_pXest)
    return new{Dc,Dp,A,B,C,D,E}(parts, dmodel, coarse_model,ptr_pXest_connectivity, ptr_pXest)
  end
end

function OctreeDistributedDiscreteModel(
  parts,
  dmodel::GridapDistributed.AbstractDistributedDiscreteModel{Dc,Dp},
  coarse_model,
  ptr_pXest_connectivity,
  ptr_pXest) where {Dc,Dp}

  return OctreeDistributedDiscreteModel(Dc,Dp,parts, dmodel, coarse_model,ptr_pXest_connectivity, ptr_pXest)
end

"""
See P4est_wrapper.jl/src/bindings/sc_common.jl for possible/valid
argument values for the p4est_verbosity_level parameter
"""
function OctreeDistributedDiscreteModel(parts::MPIData{<:Integer},
                                        coarse_model::DiscreteModel{Dc,Dp},
                                        num_uniform_refinements;
                                        p4est_verbosity_level=P4est_wrapper.SC_LP_DEFAULT) where {Dc,Dp}
  comm = parts.comm
  if i_am_in(parts.comm)
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

    return OctreeDistributedDiscreteModel(Dc,Dp,parts,dmodel,coarse_model,ptr_pXest_connectivity,ptr_pXest)
  else
    ptr_pXest_connectivity = GridapP4est.setup_pXest_connectivity(coarse_model)
    return OctreeDistributedDiscreteModel(Dc,Dp,parts,nothing,coarse_model,ptr_pXest_connectivity,nothing)
  end
end

function OctreeDistributedDiscreteModel(
    parts::MPIData{<:Integer},
    coarse_model::DiscreteModel{Dc,Dp};
    p4est_verbosity_level=P4est_wrapper.SC_LP_DEFAULT) where {Dc,Dp}
  OctreeDistributedDiscreteModel(parts,coarse_model,0; p4est_verbosity_level=p4est_verbosity_level)
end

function octree_distributed_discrete_model_free!(model::OctreeDistributedDiscreteModel{Dc}) where Dc
  if i_am_in(model.parts.comm)
    pXest_destroy(Val{Dc},model.ptr_pXest)
  end
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

GridapDistributed.local_views(model::OctreeDistributedDiscreteModel) = GridapDistributed.local_views(model.dmodel)
GridapDistributed.get_cell_gids(model::OctreeDistributedDiscreteModel) = GridapDistributed.get_cell_gids(model.dmodel)
GridapDistributed.get_face_gids(model::OctreeDistributedDiscreteModel,dim::Integer) = GridapDistributed.get_face_gids(model.dmodel,dim)
GridapDistributed.generate_gids(model::OctreeDistributedDiscreteModel,spaces) = GridapDistributed.generate_gids(model.dmodel,spaces)

###################################################################
# Private methods 

function pXest_copy(::Type{Val{Dc}}, ptr_pXest) where Dc
  if (Dc==2)
    p4est_copy(ptr_pXest, Cint(0))
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


function _create_void_octree_model(model::OctreeDistributedDiscreteModel{Dc,Dp},parts) where {Dc,Dp}
  OctreeDistributedDiscreteModel(Dc,Dp,parts,nothing,model.coarse_model,model.ptr_pXest_connectivity,nothing)
end


function _compute_fine_to_coarse_model_glue(
         cparts,
         cmodel::Union{Nothing,GridapDistributed.DistributedDiscreteModel{Dc}},
         fmodel::GridapDistributed.DistributedDiscreteModel{Dc}) where Dc

  # Fill data for owned (from coarse cells owned by the processor)
  fgids    = get_cell_gids(fmodel)
  f1,f2,f3, cgids_snd, cgids_rcv = map_parts(fmodel.models,fgids.partition) do fmodel, fpartition
    if (!(GridapP4est.i_am_in(cparts)))
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
    (fcell_to_child_id != nothing) ? fcell_to_child_id : Int[]
  end
  exchange!(dfcell_to_child_id, fgids.exchanger)
  tout = PArrays.async_exchange!(cgids_rcv,
                         cgids_snd,
                         fgids.exchanger.parts_rcv,
                         fgids.exchanger.parts_snd)
  map_parts(schedule,tout)
  map_parts(wait,tout)

  map_parts(f1,cgids_rcv) do fine_to_coarse_faces_map, cgids_rcv
    if (GridapP4est.i_am_in(cparts))
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
    if (!(GridapP4est.i_am_in(cparts)))
      nothing
    else
      reffe          = LagrangianRefFE(Float64,QUAD,1)
      ref_cell_map   = Gridap.Refinement.get_f2c_reference_cell_map(reffe,2)

      RefinementGlue(fine_to_coarse_faces_map,fcell_to_child_id,ref_cell_map)
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

function Gridap.Refinement.refine(model::OctreeDistributedDiscreteModel{Dc,Dp}, parts=nothing) where {Dc,Dp}
   comm = model.parts.comm
   if (i_am_in(comm))
     # Copy and refine input p4est
     ptr_new_pXest = pXest_copy(Val{Dc}, model.ptr_pXest)
     pXest_refine!(Val{Dc}, ptr_new_pXest)
   else
     ptr_new_pXest=nothing
   end

   comm = (parts == nothing) ? model.parts.comm : parts.comm
   if (i_am_in(comm))
      if (parts != nothing)
        aux = ptr_new_pXest
        ptr_new_pXest = _p4est_to_new_comm(ptr_new_pXest,
                                           model.ptr_pXest_connectivity,
                                           model.parts.comm,
                                           parts.comm)
        if (i_am_in(model.parts.comm))
          pXest_destroy(Val{Dc},aux)
        end
      end

      # Extract ghost and lnodes
      ptr_pXest_ghost  = setup_pXest_ghost(Val{Dc}, ptr_new_pXest)
      ptr_pXest_lnodes = setup_pXest_lnodes(Val{Dc}, ptr_new_pXest, ptr_pXest_ghost)

      # Build fine-grid mesh
      parts  = (parts == nothing) ? model.parts : parts
      fmodel = setup_distributed_discrete_model(Val{Dc},
                                              parts,
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
                                     parts,
                                     fmodel,
                                     model.coarse_model,
                                     model.ptr_pXest_connectivity,
                                     ptr_new_pXest)
      return ref_model, dglue
   else
    parts = (parts == nothing) ? model.parts : parts
    return _create_void_octree_model(model,parts), nothing
   end
end


# We have a p4est distributed among P processors. This function
# instantiates the same among Q processors.
function _p4est_to_new_comm(ptr_pXest, ptr_pXest_conn, old_comm, new_comm)
  if (GridapP4est.i_am_in(new_comm))
    new_comm_num_parts    = GridapP4est.num_parts(new_comm)
    global_first_quadrant = Vector{P4est_wrapper.p4est_gloidx_t}(undef,new_comm_num_parts+1)

    pXest_conn = ptr_pXest_conn[]
    pertree = Vector{P4est_wrapper.p4est_gloidx_t}(undef,pXest_conn.num_trees+1)
    if (GridapP4est.i_am_in(old_comm))
      pXest          = ptr_pXest[]
      old_comm_num_parts = GridapP4est.num_parts(old_comm)
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
      GridapP4est.p4est_comm_count_pertree(ptr_pXest,pertree)
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

function redistribute(model::OctreeDistributedDiscreteModel{Dc,Dp}) where {Dc,Dp}
  parts = model.parts
  if (i_am_in(parts.comm))
    ptr_pXest_old = model.ptr_pXest
    ptr_pXest     = pXest_copy(Val{Dc}, model.ptr_pXest)
    p4est_partition(ptr_pXest, 0, C_NULL)

    # Compute RedistributeGlue
    parts_snd, lids_snd, old2new = _p4est_compute_migration_control_data(Val{Dc},ptr_pXest_old,ptr_pXest)
    parts_rcv, lids_rcv, new2old = _p4est_compute_migration_control_data(Val{Dc},ptr_pXest,ptr_pXest_old)
    lids_rcv, parts_rcv = map_parts(parts) do _
      lids_rcv, parts_rcv
    end
    lids_snd, parts_snd = map_parts(parts) do _
      lids_snd, parts_snd
    end
    old2new, new2old = map_parts(parts) do _
      old2new,new2old
    end
    glue = GridapDistributed.RedistributeGlue(parts_rcv,parts_snd,lids_rcv,lids_snd,old2new,new2old)

    # Extract ghost and lnodes
    ptr_pXest_ghost  = setup_pXest_ghost(Val{Dc}, ptr_pXest)
    ptr_pXest_lnodes = setup_pXest_lnodes(Val{Dc}, ptr_pXest, ptr_pXest_ghost)

    # Build fine-grid mesh
    fmodel = setup_distributed_discrete_model(Val{Dc},
                                            model.parts,
                                            model.coarse_model,
                                            model.ptr_pXest_connectivity,
                                            ptr_pXest,
                                            ptr_pXest_ghost,
                                            ptr_pXest_lnodes)

    pXest_lnodes_destroy(Val{Dc},ptr_pXest_lnodes)
    pXest_ghost_destroy(Val{Dc},ptr_pXest_ghost)

    red_model = OctreeDistributedDiscreteModel(Dc,Dp,
                                    parts,
                                    fmodel,
                                    model.coarse_model,
                                    model.ptr_pXest_connectivity,
                                    ptr_pXest)
    return red_model, glue
  else
    return _create_void_octree_model(model,model.parts), nothing
  end
end
