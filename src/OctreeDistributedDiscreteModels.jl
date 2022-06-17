struct OctreeDistributedDiscreteModel{Dc,Dp,A,B,C,D,E} <: GridapType
  parts                  :: A
  dmodel                 :: B
  coarse_model           :: C
  ptr_pXest_connectivity :: D
  ptr_pXest              :: E
end

function OctreeDistributedDiscreteModel(parts::MPIData{<:Integer},
                              coarse_model::DiscreteModel{Dc,Dp};
                              p4est_verbosity_level=P4est_wrapper.SC_LP_DEFAULT) where {Dc,Dp}
    comm = parts.comm
    if i_am_in(parts.comm)
      sc_init(parts.comm, Cint(true), Cint(true), C_NULL, p4est_verbosity_level)
      p4est_init(C_NULL, p4est_verbosity_level)

      ptr_pXest_connectivity,
        ptr_pXest,
          ptr_pXest_ghost,
            ptr_pXest_lnodes = setup_ptr_pXest_objects(Val{Dc},
                                                        comm,
                                                        coarse_model,
                                                        0)
      dmodel=setup_distributed_discrete_model(Val{Dc},
                                              parts,
                                              coarse_model,
                                              ptr_pXest_connectivity,
                                              ptr_pXest,
                                              ptr_pXest_ghost,
                                              ptr_pXest_lnodes)
      if (Dc==2)
        # Destroy lnodes
        p4est_lnodes_destroy(ptr_pXest_lnodes)
        # Destroy ghost
        p4est_ghost_destroy(ptr_pXest_ghost)
      else
        # Destroy lnodes
        p8est_lnodes_destroy(ptr_pXest_lnodes)
        # Destroy ghost
        p8est_ghost_destroy(ptr_pXest_ghost)
      end
      A=typeof(parts)
      B=typeof(dmodel)
      C=typeof(coarse_model)
      D=typeof(ptr_pXest_connectivity)
      E=typeof(ptr_pXest)
      OctreeDistributedDiscreteModel{Dc,Dp,A,B,C,D,E}(parts,
                                                      dmodel,
                                                      coarse_model,
                                                      ptr_pXest_connectivity,
                                                      ptr_pXest)
    else
      ptr_pXest_connectivity = GridapP4est.setup_pXest_connectivity(coarse_model)
      A=typeof(parts)
      B=typeof(nothing)
      C=typeof(coarse_model)
      D=typeof(ptr_pXest_connectivity)
      E=typeof(nothing)
      OctreeDistributedDiscreteModel{Dc,Dp,A,B,C,D,E}(parts,
                                                      nothing,
                                                      coarse_model,
                                                      ptr_pXest_connectivity,
                                                      nothing)
    end
  end

function octree_distributed_discrete_model_free(model::OctreeDistributedDiscreteModel{Dc}) where Dc
  if (Dc==2)
    p4est_destroy(model.ptr_pXest)
    # p4est_connectivity_destroy(model.ptr_pXest_connectivity)
  else
    p8est_destroy(model.ptr_pXest)
    # p8est_connectivity_destroy(model.ptr_pXest_connectivity)
  end
  # sc_finalize()
end



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

struct FineToCoarseModelGlue{A,B,C}
  fine_to_coarse_faces_map::A
  fine_to_coarse_faces_dim::B
  fcell_to_child_id::C
end


function _create_void_octree_model(model::OctreeDistributedDiscreteModel{Dc,Dp}) where {Dc,Dp}
  A=typeof(model.parts)
  B=typeof(nothing)
  C=typeof(model.coarse_model)
  D=typeof(model.ptr_pXest_connectivity)
  E=typeof(nothing)
  OctreeDistributedDiscreteModel{Dc,Dp,A,B,C,D,E}(model.parts,
                                                  nothing,
                                                  model.coarse_model,
                                                  model.ptr_pXest_connectivity,
                                                  nothing)
end

function refine(model::OctreeDistributedDiscreteModel{Dc,Dp}) where {Dc,Dp}
   comm = model.parts.comm
   if (i_am_in(comm))
      # Copy and refine input p4est
      ptr_new_pXest=pXest_copy(Val{Dc}, model.ptr_pXest)
      pXest_refine!(Val{Dc}, ptr_new_pXest)

      # Extract ghost and lnodes
      ptr_pXest_ghost=setup_pXest_ghost(Val{Dc}, ptr_new_pXest)
      ptr_pXest_lnodes=setup_pXest_lnodes(Val{Dc}, ptr_new_pXest, ptr_pXest_ghost)

      # Build fine-grid mesh
      parts = get_part_ids(model.dmodel.models)
      fmodel=setup_distributed_discrete_model(Val{Dc},
                                              parts,
                                              model.coarse_model,
                                              model.ptr_pXest_connectivity,
                                              ptr_new_pXest,
                                              ptr_pXest_ghost,
                                              ptr_pXest_lnodes)

      if (Dc==2)
        # Destroy lnodes
        p4est_lnodes_destroy(ptr_pXest_lnodes)
        # Destroy ghost
        p4est_ghost_destroy(ptr_pXest_ghost)
      else
        # Destroy lnodes
        p8est_lnodes_destroy(ptr_pXest_lnodes)
        # Destroy ghost
        p8est_ghost_destroy(ptr_pXest_ghost)
        end

      dglue=map_parts(model.dmodel.models,fmodel.models) do cmodel, fmodel
        fine_to_coarse_faces_map=Vector{Vector{Int}}(undef,Dc+1)
        fine_to_coarse_faces_dim=Vector{Vector{Int}}(undef,Dc)

        num_c_cells=num_cells(cmodel)
        num_f_cells=num_cells(fmodel)

        # Allocate local vector size # local cells
        fine_to_coarse_faces_map[Dc+1]=Vector{Int}(undef,num_f_cells)
        fcell_to_child_id=Vector{Int}(undef,num_f_cells)

        # Go over all cells of coarse grid portion
        num_children=get_num_children(Val{Dc})
        c=1
        for cell=1:num_c_cells
          for child=1:num_children
            fine_to_coarse_faces_map[Dc+1][c+child-1]=cell
            fcell_to_child_id[c+child-1]=child
          end
          c=c+num_children
        end
        ctopology=Gridap.Geometry.get_grid_topology(cmodel)
        ftopology=Gridap.Geometry.get_grid_topology(fmodel)
        for d=1:Dc
          fine_to_coarse_faces_map[d] = Vector{Int}(undef,Gridap.Geometry.num_faces(ftopology,d-1))
          fine_to_coarse_faces_dim[d] = Vector{Int}(undef,Gridap.Geometry.num_faces(ftopology,d-1))
        end

        c_cell_faces=[]
        cache_c_cell_faces=[]
        f_cell_faces=[]
        cache_f_cell_faces=[]
        for d=1:Dc
          push!(c_cell_faces, Gridap.Geometry.get_faces(ctopology, Dc, d-1))
          push!(cache_c_cell_faces,array_cache(last(c_cell_faces)))
          push!(f_cell_faces, Gridap.Geometry.get_faces(ftopology, Dc, d-1))
          push!(cache_f_cell_faces,array_cache(last(f_cell_faces)))
        end
        parent_cell_faces=Vector{Vector{Int}}(undef,Dc)
        for cell=1:num_f_cells
          parent_cell=fine_to_coarse_faces_map[Dc+1][cell]
          child=fcell_to_child_id[cell]
          for d=1:Dc
            parent_cell_faces[d]=getindex!(cache_c_cell_faces[d],
                                          c_cell_faces[d],
                                          parent_cell)
          end
          for d=1:Dc
            cell_f_faces=getindex!(cache_f_cell_faces[d],
                                    f_cell_faces[d],
                                    cell)
            for (lf,f) in enumerate(cell_f_faces)
              c     = rrule_f_to_c_lid_2D[d][child][lf]
              dim_c = rrule_f_to_c_dim_2D[d][child][lf]
              if (dim_c == Dc)
                fine_to_coarse_faces_map[d][f]=parent_cell
              else
                fine_to_coarse_faces_map[d][f]=parent_cell_faces[dim_c+1][c]
              end
              fine_to_coarse_faces_dim[d][f]=dim_c
            end
          end
        end
        FineToCoarseModelGlue(fine_to_coarse_faces_map,
                              fine_to_coarse_faces_dim,
                              fcell_to_child_id)
      end
      A=typeof(model.parts)
      B=typeof(fmodel)
      C=typeof(model.coarse_model)
      D=typeof(model.ptr_pXest_connectivity)
      E=typeof(ptr_new_pXest)
      OctreeDistributedDiscreteModel{Dc,Dp,A,B,C,D,E}(model.parts,
                                      fmodel,
                                      model.coarse_model,
                                      model.ptr_pXest_connectivity,
                                      ptr_new_pXest), dglue
   else
    _create_void_octree_model(model), nothing
   end
end


# We have a p4est distributed among P processors, I want to
# instantiate the same among Q processors. How do we do that?
function _p4est_to_new_comm(ptr_pXest, ptr_pXest_conn, old_comm, new_comm)
  if (GridapP4est.i_am_in(new_comm))
    new_comm_num_parts = GridapP4est.num_parts(new_comm)
    global_first_quadrant=Vector{P4est_wrapper.p4est_gloidx_t}(undef,new_comm_num_parts+1)
    pXest_conn = ptr_pXest_conn[]
    pertree=Vector{P4est_wrapper.p4est_gloidx_t}(undef,pXest_conn.num_trees+1)
    if (GridapP4est.i_am_in(old_comm))
      pXest          = ptr_pXest[]


      old_comm_num_parts = GridapP4est.num_parts(old_comm)

      old_global_first_quadrant = unsafe_wrap(Array,
                                              pXest.global_first_quadrant,
                                              old_comm_num_parts+1)
      for i=1:length(old_global_first_quadrant)
        global_first_quadrant[i]=old_global_first_quadrant[i]
      end
      for i=length(old_global_first_quadrant)+1:length(global_first_quadrant)
        global_first_quadrant[i]=old_global_first_quadrant[end]
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
    P4est_wrapper.p4est_inflate(new_comm,
                                ptr_pXest_conn,
                                global_first_quadrant,
                                pertree,
                                quadrants,
                                C_NULL,
                                C_NULL)
  else
    nothing
  end
end

function redistribute(model::OctreeDistributedDiscreteModel{Dc,Dp}, parts) where {Dc,Dp}
  if model.parts===parts
    model
  else
    ptr_pXest = _p4est_to_new_comm(model.ptr_pXest,
                                   model.ptr_pXest_connectivity,
                                   model.parts.comm,
                                   parts.comm)
    if (GridapP4est.i_am_in(parts.comm))
      p4est_partition(ptr_pXest, 0, C_NULL)
    else
      _create_void_octree_model(model)
    end
  end
end