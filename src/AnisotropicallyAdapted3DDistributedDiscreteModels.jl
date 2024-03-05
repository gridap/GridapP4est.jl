function AnisotropicallyAdapted3DDistributedDiscreteModel(
    parts::AbstractVector{<:Integer},
    coarse_model::DiscreteModel{2,2},
    num_horizontal_uniform_refinements,
    num_vertical_uniform_refinements;
    extrusion_vector::Vector{Float64}=[0.0,0.0,1.0])

    pXest_type=P6estType()
    ptr_pXest_connectivity=setup_pXest_connectivity(pXest_type,coarse_model,extrusion_vector)
    ptr_pXest=setup_pXest(pXest_type, 
                          parts.comm, 
                          ptr_pXest_connectivity, 
                          num_horizontal_uniform_refinements,
                          num_vertical_uniform_refinements)


    ptr_pXest_ghost=setup_pXest_ghost(pXest_type,ptr_pXest)
    ptr_pXest_lnodes=setup_pXest_lnodes_nonconforming(pXest_type, ptr_pXest, ptr_pXest_ghost)

    dmodel,non_conforming_glue  = setup_non_conforming_distributed_discrete_model(pXest_type,
                                                    parts,
                                                    coarse_model,
                                                    ptr_pXest_connectivity,
                                                    ptr_pXest,
                                                    ptr_pXest_ghost,
                                                    ptr_pXest_lnodes)

    pXest_lnodes_destroy(pXest_type,ptr_pXest_lnodes)
    pXest_ghost_destroy(pXest_type,ptr_pXest_ghost)

    return OctreeDistributedDiscreteModel(3,
                                          3,
                                          parts,
                                          dmodel,
                                          non_conforming_glue,
                                          coarse_model,
                                          ptr_pXest_connectivity,
                                          ptr_pXest,
                                          pXest_type,
                                          PXestHorizontalRefinementRuleType(),
                                          true,
                                          nothing)
end

function _vertically_refine_coarsen_balance!(model::OctreeDistributedDiscreteModel{Dc,Dp}, 
                                             refinement_and_coarsening_flags::MPIArray{<:Vector}) where {Dc,Dp}

  pXest_type = model.pXest_type
  init_fn_callback_c = p6est_vertically_adapt_reset_callbacks()
  coarsen_fn_callback_c = p6est_vertically_coarsen_callbacks()
  refine_callback_c,refine_replace_callback_c = p6est_vertically_refine_callbacks()

  map(model.dmodel.models,refinement_and_coarsening_flags) do lmodel, flags
    # The length of the local flags array has to match the number of 
    # cells in the model. This includes both owned and ghost cells. 
    # Only the flags for owned cells are actually taken into account. 
    @assert num_cells(lmodel)==length(flags)
    pXest_reset_data!(pXest_type, model.ptr_pXest, Cint(sizeof(Cint)), init_fn_callback_c, pointer(flags))
  end

  # # Copy input p4est, refine and balance
  ptr_new_pXest = pXest_copy(pXest_type, model.ptr_pXest)
  p6est_vertically_refine!(ptr_new_pXest,
                           refine_callback_c,
                           refine_replace_callback_c)
  p6est_vertically_coarsen!(ptr_new_pXest, coarsen_fn_callback_c)
  pXest_balance!(pXest_type, ptr_new_pXest)
  p6est_vertically_adapt_update_flags!(model.ptr_pXest,ptr_new_pXest)
  ptr_new_pXest
end 

function vertically_adapt(model::OctreeDistributedDiscreteModel{Dc,Dp}, 
		                      refinement_and_coarsening_flags::MPIArray{<:Vector{<:Integer}};
                          parts=nothing) where {Dc,Dp}

  Gridap.Helpers.@notimplementedif parts!=nothing

  _refinement_and_coarsening_flags = map(refinement_and_coarsening_flags) do flags
    convert(Vector{Cint},flags)
  end 
  
  ptr_new_pXest = _vertically_refine_coarsen_balance!(model, _refinement_and_coarsening_flags)

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
  pXest_refinement_rule_type = PXestVerticalRefinementRuleType()
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

function setup_non_conforming_distributed_discrete_model(pXest_type::P6estType,
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


  nlvertices = map(non_conforming_glue) do ncglue
    ncglue.num_regular_faces[1]+ncglue.num_hanging_faces[1]
  end

  node_coordinates=generate_node_coordinates(pXest_type,
                                             gridap_cell_faces[1],
                                             nlvertices,
                                             ptr_pXest_connectivity,
                                             ptr_pXest,
                                             ptr_pXest_ghost)

  grid,topology=generate_grid_and_topology(pXest_type,
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

  face_labeling=generate_face_labeling(pXest_type,
                                       parts,
                                       cell_prange,
                                       coarse_discrete_model,
                                       topology,
                                       ptr_pXest,
                                       ptr_pXest_ghost)

  # _set_hanging_labels!(face_labeling,non_conforming_glue)

  discretemodel=map(grid,topology,face_labeling) do grid, topology, face_labeling
    Gridap.Geometry.UnstructuredDiscreteModel(grid,topology,face_labeling)
  end
  GridapDistributed.DistributedDiscreteModel(discretemodel,cell_prange), non_conforming_glue
end


# Extrude entity Ids:
# 1. Extrude corners to edges and corners
const corner_to_extruded_corner=[2,4,6,8]
const corner_to_extruded_edge=[1,2,3,4] 
# 2. Extrude edges to faces and edges
const edge_to_extruded_edge=[6,8,10,12]
const edge_to_extruded_face=[1,2,3,4] 
# 3. Extrude interior to upper face
const cell_to_extruded_cell=[1]
const cell_to_extruded_face=[6]

function get_bottom_to_extruded_face_lids(Db,De)
  if (Db==0 && De==0)
    corner_to_extruded_corner
  elseif (Db==0 && De==1)
    corner_to_extruded_edge
  elseif (Db==1 && De==1)
    edge_to_extruded_edge
  elseif (Db==1 && De==2)
    edge_to_extruded_face
  elseif (Db==2 && De==2)
    cell_to_extruded_face
  elseif (Db==2 && De==3)
    cell_to_extruded_cell
  else 
    @assert false 
  end
end

const corner_to_corner = [1,3,5,7]
const edge_to_edge     = [5,7,9,11]
const cell_to_face     = [5]
function get_bottom_to_overlapping_face_lids(D)
  if (D==0)
    corner_to_corner
  elseif (D==1)
    edge_to_edge
  elseif (D==2)
    cell_to_face
  else 
    @assert false 
  end
end

function transfer_entity_ids!(target_entity_ids, 
                              scell_faces,
                              tcell_faces,
                              source_to_target_lid, 
                              source_entity_ids,
                              offset)
  for (blid,bgid) in enumerate(scell_faces)
    sentity_id=source_entity_ids[bgid]
    egid=tcell_faces[source_to_target_lid[blid]]
    target_entity_ids[egid]=sentity_id+offset
  end 
end

function generate_face_labeling(pXest_type::P6estType,
                                parts,
                                cell_prange,
                                coarse_discrete_model::DiscreteModel{2,2},
                                topology,
                                ptr_pXest,
                                ptr_pXest_ghost)

  Dc=3                              

  pXest       = ptr_pXest[]
  pXest_ghost = ptr_pXest_ghost[]

  p4est_type=P4estType()
  ptr_p4est_connectivity = pXest.connectivity[].conn4
  ptr_p4est = pXest.columns
  ptr_p4est_ghost = setup_pXest_ghost(p4est_type,ptr_p4est)
  ptr_p4est_lnodes = setup_pXest_lnodes_nonconforming(p4est_type, ptr_p4est, ptr_p4est_ghost)

  bottom_boundary_model,_=setup_non_conforming_distributed_discrete_model(p4est_type,
                                                                        parts,
                                                                        coarse_discrete_model,
                                                                        ptr_p4est_connectivity,
                                                                        ptr_p4est,
                                                                        ptr_p4est_ghost,
                                                                        ptr_p4est_lnodes)

  pXest_ghost_destroy(p4est_type, ptr_p4est_ghost)
  pXest_lnodes_destroy(p4est_type, ptr_p4est_lnodes)

  bottom_boundary_model_topology = map(local_views(bottom_boundary_model)) do model 
    Gridap.Geometry.get_grid_topology(model)
  end
  
  bottom_boundary_model_labeling = map(local_views(bottom_boundary_model)) do model 
    Gridap.Geometry.get_face_labeling(model)
  end

  coarse_face_labeling = get_face_labeling(coarse_discrete_model)
  max_entity_id = _compute_max_entity_id(coarse_face_labeling)

  offset_intermediate_entities = max_entity_id
  offset_top_entities = 2*max_entity_id

  faces_to_entity=map(topology, 
                     bottom_boundary_model_topology, 
                     bottom_boundary_model_labeling) do topology, btopology, blabeling

    bcell_faces  = [Gridap.Geometry.get_faces(btopology,2,i) for i in 0:1]
    bcell_caches = [array_cache(bcell_faces[i]) for i=1:length(bcell_faces)]
    cell_faces   = [Gridap.Geometry.get_faces(topology,Dc,i) for i in 0:2]
    cell_caches  = [array_cache(cell_faces[i]) for i=1:length(cell_faces)]
    bface_entity = [Gridap.Geometry.get_face_entity(blabeling,i) for i=0:2]

    num_vertices=Gridap.Geometry.num_faces(topology,0)
    vertex_to_entity=zeros(Int,num_vertices)
    
    num_edgets=Gridap.Geometry.num_faces(topology,1)
    edget_to_entity=zeros(Int,num_edgets)
    
    num_faces=Gridap.Geometry.num_faces(topology,Dc-1)
    facet_to_entity=zeros(Int,num_faces)
    
    num_cells=Gridap.Geometry.num_faces(topology,Dc)
    cell_to_entity=zeros(Int,num_cells)

    eface_entity = [vertex_to_entity,edget_to_entity,facet_to_entity,cell_to_entity]

    gcell=1
    qcell=1

    num_trees = Cint(pXest.columns[].connectivity[].num_trees)

    # Go over trees 
    for itree = 0:num_trees-1
      tree = pXest_tree_array_index(pXest_type,pXest,itree)[]
      num_quads = Cint(tree.quadrants.elem_count)
      
      # Go over columns of current tree 
      for iquad=0:num_quads-1

        current_bcell_faces = [getindex!(bcell_caches[i],bcell_faces[i],qcell) for i=1:length(bcell_caches)]
        current_cell_faces  = [getindex!(cell_caches[i],cell_faces[i],gcell) for i=1:length(cell_caches)]
  
        # Transfer entity IDs from the bottom 2D cell to the 
        # bottom of first 3D cell in the column
        for d=0:1 
          b2oflids=get_bottom_to_overlapping_face_lids(d)
          transfer_entity_ids!(eface_entity[d+1], 
                               current_bcell_faces[d+1],
                               current_cell_faces[d+1],
                               b2oflids,
                               bface_entity[d+1],
                               0)
        end
        b2oflids=get_bottom_to_overlapping_face_lids(2)
        transfer_entity_ids!(eface_entity[3], 
                             [qcell],
                             current_cell_faces[3],
                             b2oflids,
                             bface_entity[2],
                             0)

        q = pXest_quadrant_array_index(pXest_type,tree,iquad)
        f,l=P6EST_COLUMN_GET_RANGE(q[])

        # Loop over layers within current column
        for layer=f:l-1
          current_cell_faces  = [getindex!(cell_caches[i],cell_faces[i],gcell) for i=1:length(cell_caches)]


          # Transfer entity IDs from the bottom of 3D cell to the rest of the cell by extrusion
          for db=0:2
            for de=db:db+1
              b2oflids=get_bottom_to_extruded_face_lids(db,de)
              if (layer==l-1 && db==de )
                offset=offset_top_entities
              else
                offset=offset_intermediate_entities
              end
              if (db==2 && de==3)
                transfer_entity_ids!(eface_entity[de+1], 
                           [qcell],
                           [gcell],
                           b2oflids,
                           bface_entity[db+1],
                           offset)
              elseif (db==2)
                transfer_entity_ids!(eface_entity[de+1], 
                           [qcell],
                           current_cell_faces[de+1],
                           b2oflids,
                           bface_entity[db+1],
                           offset)
              else
                @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))]  qcell=$(qcell) gcell=$(gcell) de=$(de) db=$(db) current_bcell_faces[db+1]=$(current_bcell_faces[db+1]) current_cell_faces[de+1]=$(current_cell_faces[de+1]) b2oflids=$(b2oflids) bface_entity[db+1]=$(bface_entity[db+1]) eface_entity[de+1]=$(eface_entity[de+1])"
                transfer_entity_ids!(eface_entity[de+1], 
                                     current_bcell_faces[db+1],
                                     current_cell_faces[de+1],
                                     b2oflids,
                                     bface_entity[db+1],
                                     offset)
                @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))] eface_entity[de+1]=$(eface_entity[de+1])"
              end
            end 
          end
          gcell=gcell+1
        end
        qcell=qcell+1
      end
    end 
     vertex_to_entity, edget_to_entity, facet_to_entity, cell_to_entity
 end

 vertex_to_entity  = map(x->x[1]   , faces_to_entity)
 edget_to_entity   = map(x->x[2]   , faces_to_entity)
 facet_to_entity   = map(x->x[Dc]  , faces_to_entity)
 cell_to_entity    = map(x->x[Dc+1], faces_to_entity)

 polytope = HEX



 update_face_to_entity_with_ghost_data!(vertex_to_entity,
                                        cell_prange,
                                        num_faces(polytope,0),
                                        cell_to_faces(topology,Dc,0))


 update_face_to_entity_with_ghost_data!(edget_to_entity,
                                        cell_prange,
                                        num_faces(polytope,1),
                                        cell_to_faces(topology,Dc,1))
 
 update_face_to_entity_with_ghost_data!(facet_to_entity,
                                        cell_prange,
                                        num_faces(polytope,Dc-1),
                                        cell_to_faces(topology,Dc,Dc-1))

 
 update_face_to_entity_with_ghost_data!(cell_to_entity,
                                        cell_prange,
                                        num_faces(polytope,Dc),
                                        cell_to_faces(topology,Dc,Dc))

  
faces_to_entity=[vertex_to_entity,edget_to_entity,facet_to_entity,cell_to_entity]
 
face_labeling =
  map(bottom_boundary_model_labeling,faces_to_entity...) do blabeling, faces_to_entity...
    d_to_dface_to_entity     = Vector{Vector{Int}}(undef,Dc+1)
    d_to_dface_to_entity[1]  = faces_to_entity[1]
    d_to_dface_to_entity[2]  = faces_to_entity[2]
    
    d_to_dface_to_entity[Dc]   = faces_to_entity[Dc]
    d_to_dface_to_entity[Dc+1] = faces_to_entity[Dc+1]

    bottom_interior_tag=findfirst(x->x=="interior",coarse_discrete_model.face_labeling.tag_to_name)
    @assert bottom_interior_tag != nothing
    bottom_interior_entities = coarse_discrete_model.face_labeling.tag_to_entities[bottom_interior_tag]

    intermediate_interior_entities = [e+offset_intermediate_entities for e in bottom_interior_entities]
    
    boundary_intermediate_entities_set=setdiff(
                                        Set(collect(Int32(offset_intermediate_entities+1):Int32(offset_top_entities))),
                                        Set(intermediate_interior_entities))

    boundary_intermediate_entities=[e for e in boundary_intermediate_entities_set]

    intermediate_interior_entities_set=setdiff(Set(collect(offset_intermediate_entities+1:offset_top_entities)),
                                        boundary_intermediate_entities_set)

    # Tags: bottom-boundary (1), intermediate boundary (2), top boundary (3), interior (4)  
    tag_to_entities = Vector{Vector{Int32}}(undef,4)
    tag_to_entities[1] = collect(1:max_entity_id)
    
    # All the intermediate entities except those which have been extruded from the interior of the bottom model
    tag_to_entities[2] = boundary_intermediate_entities  
    
    tag_to_entities[3] = collect(offset_top_entities+1:3*max_entity_id)
    
    # All the intermediate entities except boundary_intermediate_entities
    tag_to_entities[4] = [e for e in intermediate_interior_entities_set]
    tag_to_name = ["bottom_boundary", "intermediate_boundary", "top_boundary", "interior"] 

    face_labeling=Gridap.Geometry.FaceLabeling(d_to_dface_to_entity,
                                               tag_to_entities,
                                               tag_to_name)

    add_tag_from_tags!(face_labeling, 
        "boundary", 
        ["bottom_boundary", "intermediate_boundary", "top_boundary"])
    
    face_labeling
  end

  # map(partition(cell_prange)) do indices 
  #   println("[$(MPI.Comm_rank(MPI.COMM_WORLD))] l2g=$(local_to_global(indices))")
  #   println("[$(MPI.Comm_rank(MPI.COMM_WORLD))] l2o=$(local_to_own(indices))")
  #   println("[$(MPI.Comm_rank(MPI.COMM_WORLD))] l2p=$(local_to_owner(indices))")

    
  #   cache=PartitionedArrays.assembly_cache(indices)
    
  #   println("[$(MPI.Comm_rank(MPI.COMM_WORLD))] ns=$(cache.neighbors_snd[])")
  #   println("[$(MPI.Comm_rank(MPI.COMM_WORLD))] nr=$(cache.neighbors_rcv[])")
  #   println("[$(MPI.Comm_rank(MPI.COMM_WORLD))] li_snd=$(cache.local_indices_snd[])")
  #   println("[$(MPI.Comm_rank(MPI.COMM_WORLD))] li_rcv=$(cache.local_indices_rcv[])")
  # end

  face_labeling
end

function num_locally_owned_columns(octree_model)
  @assert octree_model.pXest_type==P6estType()
  map(octree_model.parts) do _
    pXest=octree_model.ptr_pXest[]
    num_cols = 0
    num_trees = Cint(pXest.columns[].connectivity[].num_trees)
    for itree = 0:num_trees-1
      tree = pXest_tree_array_index(octree_model.pXest_type,pXest,itree)[]
      num_cols += tree.quadrants.elem_count 
    end
    num_cols
  end
end 


function _horizontally_refine_coarsen_balance!(model::OctreeDistributedDiscreteModel{Dc,Dp}, 
                                               refinement_and_coarsening_flags::MPIArray{<:Vector}) where {Dc,Dp}

  pXest_type = model.pXest_type
  init_fn_callback_c = p6est_horizontally_adapt_reset_callbacks()
  coarsen_fn_callback_c = p6est_horizontally_coarsen_callbacks()
  refine_callback_c,refine_replace_callback_c = p6est_horizontally_refine_callbacks()

  num_cols = num_locally_owned_columns(model)

  map(refinement_and_coarsening_flags,num_cols) do flags, num_cols
    # The length of the local flags array has to match the number of locally owned columns in the model 
    @assert num_cols==length(flags)
    pXest_reset_data!(pXest_type, model.ptr_pXest, Cint(sizeof(Cint)), init_fn_callback_c, pointer(flags))
  end


  # # Copy input p4est, refine and balance
  ptr_new_pXest = pXest_copy(pXest_type, model.ptr_pXest)

  p6est_horizontally_refine!(ptr_new_pXest,
                             refine_callback_c,
                             refine_replace_callback_c)

  p6est_horizontally_coarsen!(ptr_new_pXest, coarsen_fn_callback_c)
  
  pXest_balance!(pXest_type, ptr_new_pXest)

  p6est_horizontally_adapt_update_flags!(model.ptr_pXest,ptr_new_pXest)
  
  ptr_new_pXest
end 


function horizontally_adapt(model::OctreeDistributedDiscreteModel{Dc,Dp}, 
		                        refinement_and_coarsening_flags::MPIArray{<:Vector{<:Integer}};
                            parts=nothing) where {Dc,Dp}

  Gridap.Helpers.@notimplementedif parts!=nothing

  _refinement_and_coarsening_flags = map(refinement_and_coarsening_flags) do flags
    convert(Vector{Cint},flags)
  end 
  
  ptr_new_pXest = _horizontally_refine_coarsen_balance!(model, _refinement_and_coarsening_flags)

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
  pXest_refinement_rule_type = PXestHorizontalRefinementRuleType()

  extruded_ref_coarsen_flags=map(partition(get_cell_gids(model)),refinement_and_coarsening_flags) do indices, flags
     similar(flags, length(local_to_global(indices)))
  end  

  _extrude_refinement_and_coarsening_flags!(extruded_ref_coarsen_flags,
                                            refinement_and_coarsening_flags,
                                            model.ptr_pXest,
                                            ptr_new_pXest)

  stride = pXest_stride_among_children(model.pXest_type,pXest_refinement_rule_type,model.ptr_pXest)

  adaptivity_glue = _compute_fine_to_coarse_model_glue(model.pXest_type,
                                                       pXest_refinement_rule_type,
                                                       model.parts,
                                                       model.dmodel,
                                                       fmodel,
                                                       extruded_ref_coarsen_flags,
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

function _extrude_refinement_and_coarsening_flags!(
         extruded_flags::MPIArray{<:Vector{<:Integer}},
         flags::MPIArray{<:Vector{<:Integer}},
         ptr_pXest_old,
         ptr_pXest_new)

  pXest_old  = ptr_pXest_old[]
  pXest_new  = ptr_pXest_new[]
  pXest_type = P6estType()
   
  num_trees = Cint(pXest_old.columns[].connectivity[].num_trees)
  @assert num_trees == Cint(pXest_new.columns[].connectivity[].num_trees)

  map(flags,extruded_flags) do flags,extruded_flags
    current_old_quad=1
    current_cell_old=1

    # Go over trees 
    for itree=0:num_trees-1
      tree = pXest_tree_array_index(pXest_type,pXest_old,itree)[]
      num_quads = Cint(tree.quadrants.elem_count)  
      iquad=0
      # Go over columns of current tree 
      while iquad<num_quads
        q = pXest_quadrant_array_index(pXest_type,tree,iquad)
        f,l=P6EST_COLUMN_GET_RANGE(q[])
        num_layers=l-f
        if (flags[current_old_quad]==nothing_flag)         
          for j=1:num_layers
            extruded_flags[current_cell_old] = nothing_flag 
            current_cell_old += 1
          end
          iquad+=1
          current_old_quad+=1
        elseif (flags[current_old_quad]==refine_flag)
          for j=1:num_layers
            extruded_flags[current_cell_old]=refine_flag
            current_cell_old += 1
          end
          iquad+=1
          current_old_quad+=1
        else 
          @assert flags[current_old_quad  ]==coarsen_flag
          @assert flags[current_old_quad+1]==coarsen_flag
          @assert flags[current_old_quad+2]==coarsen_flag
          @assert flags[current_old_quad+3]==coarsen_flag
          for j=1:num_layers*4
            extruded_flags[current_cell_old]=coarsen_flag
            current_cell_old += 1
          end
          iquad+=4
          current_old_quad+=4
        end
      end 
    end
  end
end