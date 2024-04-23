const P4EST_2_GRIDAP_FACET_2D  = [ 3, 4, 1, 2 ]
const GRIDAP_2_P4EST_FACET_2D  = [ 3, 4, 1, 2 ]

const P4EST_2_GRIDAP_FACET_3D  = [ 5, 6, 3, 4, 1, 2 ]
const GRIDAP_2_P4EST_FACET_3D  = [ 5, 6, 3, 4, 1, 2 ]

P4est_wrapper.quadrant_data(x::Clong) = reinterpret(P4est_wrapper.quadrant_data, x)

function p4est_get_quadrant_vertex_coordinates(connectivity::Ptr{p4est_connectivity_t},
                                               treeid::p4est_topidx_t,
                                               x::p4est_qcoord_t,
                                               y::p4est_qcoord_t,
                                               level::Int8,
                                               corner::Cint,
                                               vxy::Ptr{Cdouble})

  myself = Ref{p4est_quadrant_t}(
    p4est_quadrant_t(x,y,level,Int8(0),Int16(0),P4est_wrapper.quadrant_data(Clong(0)))
  )
  neighbour = Ref{p4est_quadrant_t}(myself[])
  if corner == 1
      p4est_quadrant_face_neighbor(myself,corner,neighbour)
  elseif corner == 2
      p4est_quadrant_face_neighbor(myself,corner+1,neighbour)
  elseif corner == 3
      p4est_quadrant_corner_neighbor(myself,corner,neighbour)
  end
  # Extract numerical coordinates of lower_left
  # corner of my corner neighbour
  p4est_qcoord_to_vertex(connectivity,
                          treeid,
                          neighbour[].x,
                          neighbour[].y,
                          vxy)
end

function p6est_get_quadrant_vertex_coordinates(connectivity::Ptr{p6est_connectivity_t},
                                               treeid::p4est_topidx_t,
                                               x::p4est_qcoord_t,
                                               y::p4est_qcoord_t,
                                               z::p4est_qcoord_t,
                                               xylevel::Int8,
                                               zlevel::Int8,
                                               corner::Cint,
                                               pvxyz::Ptr{Cdouble})
   
   function p6est_corner_to_p4est_corner(p6est_corner::Cint)
      div(p6est_corner,Cint(2))
   end 

   function zquadrant_to_zcorner(z, zlevel_quadrant, p6est_corner)
      if p6est_corner in (0,2,4,6)
        z
      elseif p6est_corner in (1,3,5,7)
        z + P4EST_QUADRANT_LEN(zlevel)
      end 
   end

   p6est_qcoord_to_vertex(connectivity,
                          treeid,
                          x, 
                          y,
                          zquadrant_to_zcorner(z,zlevel,corner), 
                          pvxyz)

   vxyz=unsafe_wrap(Array, pvxyz, 3)
   zcoord=vxyz[3]
   
   # Always sets the z-coordinate to 0.0!
   p4est_get_quadrant_vertex_coordinates(connectivity[].conn4,
                                         treeid,
                                         x,
                                         y,
                                         xylevel,
                                         p6est_corner_to_p4est_corner(corner),
                                         pvxyz)
   vxyz[3]=zcoord
end

function  p8est_get_quadrant_vertex_coordinates(connectivity::Ptr{p8est_connectivity_t},
                                                treeid::p4est_topidx_t,
                                                x::p4est_qcoord_t,
                                                y::p4est_qcoord_t,
                                                z::p4est_qcoord_t,
                                                level::Int8,
                                                corner::Cint,
                                                vxyz::Ptr{Cdouble})

  myself = Ref{p8est_quadrant_t}(
    p8est_quadrant_t(x,y,z,level,Int8(0),Int16(0),P4est_wrapper.quadrant_data(Clong(0)))
  )
  neighbour = Ref{p8est_quadrant_t}(myself[])

  if ( corner == 1 )
    p8est_quadrant_face_neighbor(myself,Cint(1),neighbour)
  elseif ( corner == 2 )
    p8est_quadrant_face_neighbor(myself,Cint(3),neighbour)
  elseif ( corner == 3 )
    p8est_quadrant_edge_neighbor(myself,Cint(11),neighbour)
  elseif ( corner == 4 )
    p8est_quadrant_face_neighbor(myself,Cint(5),neighbour)
  elseif ( corner == 5 )
    p8est_quadrant_edge_neighbor(myself,Cint(7),neighbour)
  elseif ( corner == 6 )
    p8est_quadrant_edge_neighbor(myself,Cint(3),neighbour)
  elseif ( corner == 7 )
    p8est_quadrant_corner_neighbor(myself,Cint(7),neighbour)
  end

  # Extract numerical coordinates of lower_left corner of my corner neighbour
  p8est_qcoord_to_vertex(connectivity,
                         treeid,
                         neighbour[].x,
                         neighbour[].y,
                         neighbour[].z,
                         vxyz)
end

"""
    setup_pXest_connectivity(cmodel::DiscreteModel)
  
  Returns the pXest connectivity object from a serial DiscreteModel.

  Two cases:

    - n_vertices == n_corners: Geometrical and topological information are the same.
    - n_vertices != n_corners: Geometrical and topological information are different. This 
                               is the case for periodic models.
"""
function setup_pXest_connectivity(cmodel::DiscreteModel)
  n_vertices = num_nodes(cmodel)
  n_corners = num_vertices(cmodel)

  if n_vertices == n_corners
    return setup_pXest_connectivity_from_geometry(cmodel)
  else
    return setup_pXest_connectivity_with_topology(cmodel)
  end
end

"""
    setup_pXest_connectivity_from_geometry(cmodel::DiscreteModel{Dc,Dp})
  
  Returns the pXest connectivity object from a serial DiscreteModel, assuming that geometrical
  and topological information are the same.
"""
function setup_pXest_connectivity_from_geometry(cmodel::DiscreteModel{Dc,Dp}) where {Dc,Dp}
  n_trees = num_cells(cmodel)
  n_vertices = num_nodes(cmodel)

  if (Dc==2)
    pconn = p4est_connectivity_new(
        p4est_topidx_t(n_vertices), # num_vertices
        p4est_topidx_t(n_trees),    # num_trees
        p4est_topidx_t(0),
        p4est_topidx_t(0))
  else
    @assert Dc==3
    pconn = p8est_connectivity_new(
        p4est_topidx_t(n_vertices), # num_vertices
        p4est_topidx_t(n_trees),    # num_trees
        p4est_topidx_t(0),
        p4est_topidx_t(0),
        p4est_topidx_t(0),
        p4est_topidx_t(0))
  end
  conn = pconn[]

  # Fill geometrical information, i.e `vertices` and `tree_to_vertex`.
  vertex_coordinates = Gridap.Geometry.get_node_coordinates(cmodel)
  vertices = unsafe_wrap(Array, conn.vertices, n_vertices*3)
  for (i,p) in enumerate(vertex_coordinates)
    offset = (i-1)*3
    vertices[offset+1] = Cdouble(p[1])
    vertices[offset+2] = Cdouble(p[2])
    vertices[offset+3] = (Dp==3) ? Cdouble(p[3]) : Cdouble(0.0) # Z coordinate always to 0.0 in 2D
  end

  cell_to_vertex = Gridap.Geometry.get_cell_node_ids(cmodel)
  tree_to_vertex = unsafe_wrap(Array, conn.tree_to_vertex, n_trees*(2^Dc))
  c2v_cache = Gridap.Arrays.array_cache(cell_to_vertex)
  for cell = 1:n_trees
    offset = (cell-1)*(2^Dc)
    vertices = Gridap.Arrays.getindex!(c2v_cache,cell_to_vertex,cell)
    for (i,v) in enumerate(vertices)
      tree_to_vertex[offset+i] = p4est_topidx_t(v-1)
    end
  end

  # Mockup `tree_to_tree` and `tree_to_face` to make sure we have a valid connectivity.
  PXEST_FACES = 2*Dc
  tree_to_tree = unsafe_wrap(Array, conn.tree_to_tree, conn.num_trees*PXEST_FACES )
  tree_to_face = unsafe_wrap(Array, conn.tree_to_face, conn.num_trees*PXEST_FACES )
  for tree = 1:conn.num_trees
    for face = 1:PXEST_FACES
      tree_to_tree[PXEST_FACES * (tree-1) + face] = tree-1
      tree_to_face[PXEST_FACES * (tree-1) + face] = face-1
    end
  end

  # Let p4est/p8est fill the topological information automatically.
  if (Dc==2)
    p4est_connectivity_complete(pconn)
    @assert Bool(p4est_connectivity_is_valid(pconn))
  else
    p8est_connectivity_complete(pconn)
    @assert Bool(p8est_connectivity_is_valid(pconn))
  end

  return pconn
end

"""
    setup_pXest_connectivity_with_topology(cmodel::DiscreteModel{Dc,Dp})
  
  Returns the pXest connectivity object from a serial DiscreteModel, manually completing 
  both geometrical and topological information. This is needed if geometry != topology, which 
  is the case for periodic models.
"""
function setup_pXest_connectivity_with_topology(cmodel::DiscreteModel{Dc,Dp}) where {Dc,Dp}
  
  # In periodic models, n_vertices != n_corners
  topo = Gridap.Geometry.get_grid_topology(cmodel)
  n_trees = num_cells(cmodel)
  n_vertices = num_nodes(cmodel)
  n_corners = Gridap.Geometry.num_faces(topo,0)
  n_edges = Gridap.Geometry.num_faces(topo,1) # In 2D, == n_facets
  n_ctt = length(Gridap.Geometry.get_faces(topo,0,Dc).data)
  n_ett = length(Gridap.Geometry.get_faces(topo,1,Dc).data)

  if (Dc==2)
    pconn = p4est_connectivity_new(
        p4est_topidx_t(n_vertices), # num_vertices
        p4est_topidx_t(n_trees),    # num_trees
        p4est_topidx_t(n_corners),  # num_corners
        p4est_topidx_t(n_ctt))
  else
    @assert Dc==3
    pconn = p8est_connectivity_new(
        p4est_topidx_t(n_vertices), # num_vertices
        p4est_topidx_t(n_trees),    # num_trees
        p4est_topidx_t(n_edges),    # num_edges
        p4est_topidx_t(n_ett),
        p4est_topidx_t(n_corners),  # num_corners
        p4est_topidx_t(n_ctt))
  end
  conn = pconn[]

  # A) Fill geometrical information, i.e `vertices` and `tree_to_vertex`.
  vertex_coordinates = Gridap.Geometry.get_node_coordinates(cmodel)
  vertices = unsafe_wrap(Array, conn.vertices, n_vertices*3)
  for (i,p) in enumerate(vertex_coordinates)
    offset = (i-1)*3
    vertices[offset+1] = Cdouble(p[1])
    vertices[offset+2] = Cdouble(p[2])
    vertices[offset+3] = (Dp==3) ? Cdouble(p[3]) : Cdouble(0.0) # Z coordinate always to 0.0 in 2D
  end

  cell_to_vertex = Gridap.Geometry.get_cell_node_ids(cmodel)
  tree_to_vertex = unsafe_wrap(Array, conn.tree_to_vertex, n_trees*(2^Dc))
  c2v_cache = Gridap.Arrays.array_cache(cell_to_vertex)
  for cell = 1:n_trees
    offset = (cell-1)*(2^Dc)
    vertices = Gridap.Arrays.getindex!(c2v_cache,cell_to_vertex,cell)
    for (i,v) in enumerate(vertices)
      tree_to_vertex[offset+i] = p4est_topidx_t(v-1)
    end
  end

  # B) Fill corner topological information
  # -> `tree_to_corner` : Topological equivalent to `tree_to_vertex`.
  # -> `corner_to_tree` and `corner_to_corner` : 
  #     For each corner, we have that `corner_to_X[ctt_offset[corner]:ctt_offset[corner+1]-1]`
  #     contains
  #      - the cell ids for the sourrounding cells
  #      - its local corner id for those cells
  cell_to_corner = Gridap.Geometry.get_faces(topo,Dc,0)
  tree_to_corner = unsafe_wrap(Array, conn.tree_to_corner, n_trees*(2^Dc))
  t2c_cache = Gridap.Arrays.array_cache(cell_to_corner)
  for cell = 1:n_trees
    offset = (cell-1)*(2^Dc)
    corners = Gridap.Arrays.getindex!(t2c_cache,cell_to_corner,cell)
    for (i,c) in enumerate(corners)
      tree_to_corner[offset+i] = p4est_topidx_t(c-1)
    end
  end

  corner_to_cell = Gridap.Geometry.get_faces(topo,0,Dc)
  ctt_offset = unsafe_wrap(Array, conn.ctt_offset, n_corners+1)
  corner_to_tree = unsafe_wrap(Array, conn.corner_to_tree, n_ctt)
  corner_to_corner = unsafe_wrap(Array, conn.corner_to_corner, n_ctt)
  ctt_offset[1] = 0
  c2t_cache = Gridap.Arrays.array_cache(corner_to_cell)
  for corner in 1:n_corners
    offset = ctt_offset[corner]
    cells = Gridap.Arrays.getindex!(c2t_cache,corner_to_cell,corner)
    ctt_offset[corner+1] = offset + length(cells)
    for (i,c) in enumerate(cells)
      cell_corners = getindex!(t2c_cache,cell_to_corner,c)
      corner_lid = findfirst(x->x==corner,cell_corners)
      corner_to_tree[offset+i] = p4est_topidx_t(c-1)
      corner_to_corner[offset+i] = Int8(corner_lid-1)
    end
  end
  @assert ctt_offset[n_corners+1] == n_ctt

  # C) Fill edge topological information. Only for 3D. 
  # Variables are equivalent to corner ones.
  if Dc == 3
    cell_to_edge = Gridap.Geometry.get_faces(topo,Dc,1)
    tree_to_edge = unsafe_wrap(Array, conn.tree_to_edge, n_trees*12)
    t2e_cache = Gridap.Arrays.array_cache(cell_to_edge)
    for cell = 1:n_trees
      offset = (cell-1)*12
      edges = Gridap.Arrays.getindex!(t2e_cache,cell_to_edge,cell)
      for (i,e) in enumerate(edges)
        tree_to_edge[offset+i] = p4est_topidx_t(e-1)
      end
    end

    edge_to_cell = Gridap.Geometry.get_faces(topo,1,Dc)
    ett_offset = unsafe_wrap(Array, conn.ett_offset, n_edges+1)
    edge_to_tree = unsafe_wrap(Array, conn.edge_to_tree, n_ett)
    edge_to_edge = unsafe_wrap(Array, conn.edge_to_edge, n_ett)
    ett_offset[1] = 0
    e2t_cache = Gridap.Arrays.array_cache(edge_to_cell)
    for edge in 1:n_edges
      offset = ett_offset[edge]
      cells = Gridap.Arrays.getindex!(e2t_cache,edge_to_cell,edge)
      ett_offset[edge+1] = offset + length(cells)
      for (i,c) in enumerate(cells)
        cell_edges = getindex!(t2e_cache,cell_to_edge,c)
        edge_lid = findfirst(x->x==edge,cell_edges)
        edge_to_tree[offset+i] = p4est_topidx_t(c-1)
        edge_to_edge[offset+i] = Int8(edge_lid-1)
      end
    end
    @assert ett_offset[n_edges+1] == n_ett
  end

  # D) Fill cell and face topological information, i.e `tree_to_tree` and `tree_to_face`.
  # Given a cell and a local face id, tree_to_X[cell*PXEST_FACES+face] contains: 
  #   - the nboring cell id for that face, or the cell id itself if it is a boundary face.
  #   - the local face id of that face in the neighboring cell, plus the relative orientation
  #     of the face in the neighboring cell. 
  #     See [p4est_connectivity.h](https://p4est.github.io/api/p4est-latest/structp4est__connectivity.html) 
  #     for more details.
  PXEST_FACES = 2*Dc
  face_to_cell = Gridap.Geometry.get_faces(topo,Dc-1,Dc)
  cell_to_face = Gridap.Geometry.get_faces(topo,Dc,Dc-1)
  tree_to_tree = unsafe_wrap(Array, conn.tree_to_tree, conn.num_trees*PXEST_FACES)
  tree_to_face = unsafe_wrap(Array, conn.tree_to_face, conn.num_trees*PXEST_FACES)

  GRIDAP_TO_P4EST_FACE = (Dc==2) ? GRIDAP_2_P4EST_FACET_2D : GRIDAP_2_P4EST_FACET_3D
  pXest_face_corners = (Dc==2) ? p4est_face_corners : p8est_face_corners

  f2t_cache = Gridap.Arrays.array_cache(face_to_cell)
  t2f_cache = Gridap.Arrays.array_cache(cell_to_face)
  nbor_t2c_cache = Gridap.Arrays.array_cache(cell_to_corner)
  nbor_t2f_cache = Gridap.Arrays.array_cache(cell_to_face)
  for tree = 1:conn.num_trees
    corners = Gridap.Arrays.getindex!(t2c_cache,cell_to_corner,tree)
    faces = view(Gridap.Arrays.getindex!(t2f_cache,cell_to_face,tree),GRIDAP_TO_P4EST_FACE)
    @assert length(faces) == PXEST_FACES
    for (face_lid,face) in enumerate(faces)
      face_cells = Gridap.Arrays.getindex!(f2t_cache,face_to_cell,face)
      nbor_lid = findfirst(x->x!=tree,face_cells)
      nbor = isnothing(nbor_lid) ? tree : face_cells[nbor_lid]
      tree_to_tree[PXEST_FACES * (tree-1) + face_lid] = p4est_topidx_t(nbor-1)

      nbor_corners = Gridap.Arrays.getindex!(nbor_t2c_cache,cell_to_corner,nbor)
      nbor_faces = view(Gridap.Arrays.getindex!(nbor_t2f_cache,cell_to_face,nbor),GRIDAP_TO_P4EST_FACE)
      nbor_face_lid = findfirst(x->x==face,nbor_faces)

      face_corners = corners[pXest_face_corners[face_lid,:].+1]
      nbor_face_corners = nbor_corners[pXest_face_corners[nbor_face_lid,:].+1]

      face_orientation = 0
      if face_lid < nbor_face_lid
        face_orientation = findfirst(x->x==face_corners[1],nbor_face_corners) - 1
      else
        face_orientation = findfirst(x->x==nbor_face_corners[1],face_corners) - 1
      end
      tree_to_face[PXEST_FACES * (tree-1) + face_lid] = Int8(PXEST_FACES*face_orientation + nbor_face_lid-1)
    end
  end

  # Make sure we have a valid connectivity
  if (Dc==2)
    @assert Bool(p4est_connectivity_is_valid(pconn))
  else
    @assert Bool(p8est_connectivity_is_valid(pconn))
  end

  return pconn
end

function fetch_vector_ghost_values_cache(vector_partition,partition)
  cache = PArrays.p_vector_cache(vector_partition,partition)
  map(reverse,cache)
end 

function fetch_vector_ghost_values!(vector_partition,cache)
  assemble!((a,b)->b, vector_partition, cache) 
end

function generate_cell_vertex_gids(ptr_pXest_lnodes, cell_prange)
  pXest_lnodes = ptr_pXest_lnodes[]

  nvertices = pXest_lnodes.vnodes
  element_nodes = unsafe_wrap(Array,
                              pXest_lnodes.element_nodes,
                              pXest_lnodes.num_local_elements*nvertices)

  nonlocal_nodes = unsafe_wrap(Array,
                               pXest_lnodes.nonlocal_nodes,
                               pXest_lnodes.num_local_nodes-pXest_lnodes.owned_count)

  k = 1
  cell_vertex_gids = map(partition(cell_prange)) do indices
    n = length(local_to_own(indices))
    ptrs = Vector{Int32}(undef,n+1)
    ptrs[1] = 1
    for i = 1:n
      ptrs[i+1] = ptrs[i] + nvertices
    end
    k = 1
    current = 1
    data = Vector{Int}(undef, ptrs[n+1]-1)
    for i = 1:pXest_lnodes.num_local_elements
      for j = 1:nvertices
        l = element_nodes[k+j-1]
        if (l < pXest_lnodes.owned_count)
          data[current] = pXest_lnodes.global_offset+l+1
        else
          data[current] = nonlocal_nodes[l-pXest_lnodes.owned_count+1]+1
        end
        current = current+1
      end
      k = k + nvertices
    end
    PArrays.JaggedArray(data,ptrs)
  end
  fetch_cache = fetch_vector_ghost_values_cache(cell_vertex_gids,partition(cell_prange))
  fetch_vector_ghost_values!(cell_vertex_gids,fetch_cache) |> wait
  return cell_vertex_gids
end

"""
  Generate the topological cellwise corner lids from their gids.
"""
function generate_cell_corner_lids(cell_corner_gids)
  map(cell_corner_gids) do cell_corner_gids
    g2l = Dict{Int,Int}()
    lid = 1
    data = Vector{Int}(undef,length(cell_corner_gids.data))
    for (i,gid) in enumerate(cell_corner_gids.data)
      if haskey(g2l,gid)
        data[i] = g2l[gid]
      else
        data[i] = lid
        g2l[gid] = lid
        lid = lid+1
      end
    end
    Gridap.Arrays.Table(data,cell_corner_gids.ptrs)
  end
end

"""
  Generate the geometrical cellwise vertex coordinates.
"""
function generate_cell_vertex_coordinates(pXest_type::PXestType,
                                          cell_vertex_lids,
                                          ptr_pXest_connectivity,
                                          ptr_pXest,
                                          ptr_pXest_ghost)

  Dc = num_cell_dims(pXest_type)

  PXEST_CORNERS = 2^Dc
  pXest_ghost = ptr_pXest_ghost[]
  pXest       = ptr_pXest[]

  cell_vertex_coordinates = map(cell_vertex_lids) do cell_vertex_lids
    data = Vector{Point{Dc,Float64}}(undef,length(cell_vertex_lids.data))
    current = 1
    vxy = Vector{Cdouble}(undef,Dc)
    pvxy = pointer(vxy,1)
    for itree = 1:pXest_ghost.num_trees
      tree = pXest_tree_array_index(pXest_type, pXest, itree-1)[]
      for cell=1:tree.quadrants.elem_count
        quadrant=pXest_quadrant_array_index(pXest_type,tree,cell-1)[] 
        # Loop over layers in the current column
        for l=1:pXest_num_quadrant_layers(pXest_type,quadrant)
          layer=pXest_get_layer(pXest_type, quadrant, pXest, l-1)
          coords=pXest_cell_coords(pXest_type,quadrant,layer)
          levels=pXest_get_quadrant_and_layer_levels(pXest_type,quadrant,layer)
          
          for vertex=1:PXEST_CORNERS
            pXest_get_quadrant_vertex_coordinates(pXest_type,
                                                  ptr_pXest_connectivity,
                                                  p4est_topidx_t(itree-1),
                                                  coords,
                                                  levels,
                                                  Cint(vertex-1),
                                                  pvxy)
            data[current] = Point{Dc,Float64}(vxy...)
            current=current+1
          end
        end 
      end
    end

    tree_offsets = unsafe_wrap(Array, pXest_ghost.tree_offsets, pXest_ghost.num_trees+1)
    ptr_ghost_quadrants = _unwrap_ghost_quadrants(pXest_type, pXest_ghost)

    # Go over ghost cells
    for i=1:pXest_ghost.num_trees
      for j=tree_offsets[i]:tree_offsets[i+1]-1
          quadrant = ptr_ghost_quadrants[j+1]
          coords=pXest_cell_coords(pXest_type,quadrant,0)
          levels=pXest_get_quadrant_and_layer_levels(pXest_type,quadrant,0)
          for vertex=1:PXEST_CORNERS
              pXest_get_quadrant_vertex_coordinates(pXest_type,
                                                    ptr_pXest_connectivity,
                                                    p4est_topidx_t(i-1),
                                                    coords,
                                                    levels,
                                                    Cint(vertex-1),
                                                    pvxy)
              data[current]=Point{Dc,Float64}(vxy...)
              current=current+1
          end
      end
    end
    ptrs = copy(cell_vertex_lids.ptrs)
    return Gridap.Arrays.Table(data,ptrs)
  end
  return cell_vertex_coordinates
end

"""
  From the topological cellwise corner ids and the geometrical cellwise vertex coordinates,
  generate the remaining topological and geometrical information. 
  I.e 
   - the geometrical vertex coordinates and the cell-wise vertex ids
   - the topological corner coordinates
"""
function generate_coords(
  topo_cell_ids :: Gridap.Arrays.Table{<:Ti},
  model_cell_coords :: Gridap.Arrays.Table{<:VectorValue{Dp,T}}
) where {Ti,Dp,T}
  n_corners = maximum(topo_cell_ids.data;init=0)
  
  n_vertices_max = length(model_cell_coords.data)
  model_coords = fill(VectorValue(fill(T(Inf),Dp)),n_vertices_max)
  for (vertex,coord) in zip(topo_cell_ids.data,model_cell_coords.data)
    model_coords[vertex] = min(model_coords[vertex],coord)
  end
  
  n_vertices = n_corners
  model_cell_ids = Gridap.Arrays.Table(copy(topo_cell_ids.data),copy(topo_cell_ids.ptrs))
  for (k,(vertex,coord)) in enumerate(zip(topo_cell_ids.data,model_cell_coords.data))
    if norm(coord-model_coords[vertex]) > eps(T)
      pos = findfirst(x -> norm(x-coord) < eps(T), model_coords[n_corners+1:n_vertices])
      if !isnothing(pos)
        model_cell_ids.data[k] = model_coords[n_corners+pos]
      else
        n_vertices += 1
        model_coords[n_vertices] = coord
        model_cell_ids.data[k] = n_vertices
      end
    end
  end
  @debug "function generate_coords :: n_vertices=$n_vertices, n_corners=$n_corners."

  resize!(model_coords,n_vertices)
  topo_coords = (n_vertices == n_corners) ? model_coords : model_coords[1:n_corners]
  return model_cell_ids, model_coords, topo_coords
end

"""
  Generate a new grid and topology from 
    - the cellwise corner ids (topology)
    - the cellwise vertex coordinates (geometry)
"""
function generate_grid_and_topology(pXest_type::P4P8estType,
                                    cell_corner_lids,
                                    cell_vertex_coordinates)

  Dc = num_cell_dims(pXest_type)

  map(cell_corner_lids, cell_vertex_coordinates) do cell_corner_lids, cell_vertex_coordinates

    cell_vertex_lids, vertex_coords, corner_coords = generate_coords(
      cell_corner_lids, cell_vertex_coordinates
    )

    poly  = (Dc==2) ? QUAD : HEX
    reffe = Gridap.ReferenceFEs.ReferenceFE(poly,lagrangian,Float64,1)
    cell_types = fill(1,length(cell_vertex_lids))

    grid = Gridap.Geometry.UnstructuredGrid(
      vertex_coords,cell_vertex_lids,[reffe],cell_types,Gridap.Geometry.NonOriented()
    )
    topology = Gridap.Geometry.UnstructuredGridTopology(
      corner_coords,cell_corner_lids,cell_types,[poly],Gridap.Geometry.NonOriented()
    )
    return grid, topology
  end |> tuple_of_arrays
end

const ITERATOR_RESTRICT_TO_BOUNDARY=Cint(100)
const ITERATOR_RESTRICT_TO_INTERIOR=Cint(101)

function generate_face_labeling(pXest_type::P4P8estType,
                                parts,
                                cell_prange,
                                coarse_discrete_model::DiscreteModel{Dc,Dp},
                                topology,
                                ptr_pXest,
                                ptr_pXest_ghost) where {Dc,Dp}

  pXest       = ptr_pXest[]
  pXest_ghost = ptr_pXest_ghost[]

  coarse_grid_topology = Gridap.Geometry.get_grid_topology(coarse_discrete_model)
  coarse_grid_labeling = Gridap.Geometry.get_face_labeling(coarse_discrete_model)
  coarse_cell_vertices = Gridap.Geometry.get_faces(coarse_grid_topology,Dc,0)
  if (Dc==3)
    coarse_cell_edgets = Gridap.Geometry.get_faces(coarse_grid_topology,Dc,1)
  end
  coarse_cell_facets   = Gridap.Geometry.get_faces(coarse_grid_topology,Dc,Dc-1)

  owned_trees_offset=Vector{Int}(undef,pXest_ghost.num_trees+1)
  owned_trees_offset[1]=0
  for itree=1:pXest_ghost.num_trees
    if Dc==2
      tree = p4est_tree_array_index(pXest.trees, itree-1)[]
    else
      tree = p8est_tree_array_index(pXest.trees, itree-1)[]
    end
    owned_trees_offset[itree+1]=owned_trees_offset[itree]+tree.quadrants.elem_count
  end

 faces_to_entity=map(topology) do topology
     # Iterate over corners
     num_vertices=Gridap.Geometry.num_faces(topology,0)
     vertex_to_entity=zeros(Int,num_vertices)
     cell_vertices=Gridap.Geometry.get_faces(topology,Dc,0)

     # Corner iterator callback
     function jcorner_callback(pinfo     :: Ptr{p8est_iter_corner_info_t},
                               user_data :: Ptr{Cvoid})
        info=pinfo[]
        if (Dc==2)
          sides=Ptr{p4est_iter_corner_side_t}(info.sides.array)
          CONNECT_CORNER=P4est_wrapper.P4EST_CONNECT_CORNER
        else
          sides=Ptr{p8est_iter_corner_side_t}(info.sides.array)
          CONNECT_CORNER=P4est_wrapper.P8EST_CONNECT_CORNER
        end
        nsides=info.sides.elem_count
        tree=sides[1].treeid+1
        data=sides[1]
        if data.is_ghost==1
           ref_cell=pXest.local_num_quadrants+data.quadid+1
        else
           ref_cell=owned_trees_offset[tree]+data.quadid+1
        end
        corner=sides[1].corner+1
        ref_cornergid=cell_vertices[ref_cell][corner]
        if (info.tree_boundary!=0 && info.tree_boundary==CONNECT_CORNER)
              # The current corner is also a corner of the coarse mesh
              coarse_cornergid=coarse_cell_vertices[tree][corner]
              vertex_to_entity[ref_cornergid]=
                 coarse_grid_labeling.d_to_dface_to_entity[1][coarse_cornergid]
              @debug "[GLOBAL CORNER] vertex_to_entity[$(ref_cornergid)]=$(coarse_grid_labeling.d_to_dface_to_entity[1][coarse_cornergid])"
        else
          if vertex_to_entity[ref_cornergid]==0
            # We are on the interior of a tree (if we did not touch it yet)
            vertex_to_entity[ref_cornergid]=coarse_grid_labeling.d_to_dface_to_entity[Dc+1][tree]
            @debug "[INTERIOR CORNER] vertex_to_entity[$(ref_cornergid)]=$(coarse_grid_labeling.d_to_dface_to_entity[Dc+1][tree])"
          end
        end
        nothing
     end

     #  C-callable face callback
     ccorner_callback = @cfunction($jcorner_callback,
                                   Cvoid,
                                   (Ptr{p8est_iter_corner_info_t},Ptr{Cvoid}))

     cell_edgets=Gridap.Geometry.get_faces(topology,Dc,1)
     num_edgets=Gridap.Geometry.num_faces(topology,1)
     edget_to_entity=zeros(Int,num_edgets)
    if (Dc==3)
      # Edge iterator callback
      function jedge_callback(pinfo     :: Ptr{p8est_iter_edge_info_t},
                               user_data :: Ptr{Cvoid})
        info=pinfo[]
        sides=Ptr{p8est_iter_edge_side_t}(info.sides.array)
        function process_edget(tree,edge,ref_cell,info)
           polytope=HEX
           poly_faces=Gridap.ReferenceFEs.get_faces(polytope)
           poly_edget_range=Gridap.ReferenceFEs.get_dimrange(polytope,1)
           poly_first_edget=first(poly_edget_range)
           poly_facet=poly_first_edget+edge-1
           if (info.tree_boundary!=0 && info.tree_boundary==P4est_wrapper.P8EST_CONNECT_EDGE)
             coarse_edgetgid=coarse_cell_edgets[tree][edge]
             coarse_edgetgid_entity=coarse_grid_labeling.d_to_dface_to_entity[2][coarse_edgetgid]
             # We are on the boundary of coarse mesh or inter-octree boundary
             for poly_incident_face in poly_faces[poly_facet]
               if poly_incident_face == poly_facet
                 ref_edgetgid=cell_edgets[ref_cell][edge]
                 edget_to_entity[ref_edgetgid]=coarse_edgetgid_entity
               else
                 ref_cornergid=cell_vertices[ref_cell][poly_incident_face]
                 vertex_to_entity[ref_cornergid]=coarse_edgetgid_entity
               end
             end
           else
             # We are on the interior of the domain if we did not touch the edge yet
             ref_edgetgid=cell_edgets[ref_cell][edge]
             if (edget_to_entity[ref_edgetgid]==0)
               edget_to_entity[ref_edgetgid]=coarse_grid_labeling.d_to_dface_to_entity[Dc+1][tree]
             end
           end
        end 
  
        nsides=info.sides.elem_count
        for iside=1:nsides
         edge=sides[iside].edge+1
         tree=sides[iside].treeid+1
         if (sides[iside].is_hanging == 0)
           data=sides[iside].is.full
           if data.is_ghost==1
             ref_cell=pXest.local_num_quadrants+data.quadid+1
           else
             ref_cell=owned_trees_offset[tree]+data.quadid+1
           end 
           process_edget(tree,edge,ref_cell,info)
         else 
           for i=1:length(sides[iside].is.hanging.quadid)
             quadid=sides[iside].is.hanging.quadid[i]
             if (sides[iside].is.hanging.is_ghost[i]==1)
               ref_cell=pXest.local_num_quadrants+quadid+1
             else 
               ref_cell=owned_trees_offset[tree]+quadid+1
             end 
             process_edget(tree,edge,ref_cell,info)
           end 
         end 
        end
        nothing
      end

      # C-callable edge callback
      cedge_callback = @cfunction($jedge_callback,
                                   Cvoid,
                                   (Ptr{p8est_iter_edge_info_t},Ptr{Cvoid}))
    end

    # Iterate over faces
    num_faces=Gridap.Geometry.num_faces(topology,Dc-1)
    facet_to_entity=zeros(Int,num_faces)
    cell_facets=Gridap.Geometry.get_faces(topology,Dc,Dc-1)

     # Face iterator callback
     function jface_callback(pinfo     :: Ptr{p8est_iter_face_info_t},
                             user_data :: Ptr{Cvoid})
        info=pinfo[]
        if Dc==2
          sides=Ptr{p4est_iter_face_side_t}(info.sides.array)
        else
          sides=Ptr{p8est_iter_face_side_t}(info.sides.array)
        end
        ptr_user_data=Ptr{Cint}(user_data)
        iterator_mode=unsafe_wrap(Array, ptr_user_data, 1)
        if (iterator_mode[1]==ITERATOR_RESTRICT_TO_BOUNDARY)
          # If current face is NOT in the boundary
          if (info.tree_boundary==0)
            return nothing
          end
        else
          # If current face is in the boundary 
          if (info.tree_boundary!=0)
            return nothing
          end
        end

        function process_facet(tree,face,ref_cell,info)
          if Dc==2
            gridap_facet=P4EST_2_GRIDAP_FACET_2D[face]
          else
            gridap_facet=P4EST_2_GRIDAP_FACET_3D[face]
          end

          polytope= Dc==2 ? QUAD : HEX
          poly_faces=Gridap.ReferenceFEs.get_faces(polytope)
          poly_facet_range=Gridap.ReferenceFEs.get_dimrange(polytope,Dc-1)
          poly_first_facet=first(poly_facet_range)
          poly_facet=poly_first_facet+gridap_facet-1

          if (info.tree_boundary!=0)
            # We are on the boundary of coarse mesh or inter-octree boundary
            coarse_facetgid=coarse_cell_facets[tree][gridap_facet]
            coarse_facetgid_entity=coarse_grid_labeling.d_to_dface_to_entity[Dc][coarse_facetgid]
          else
            coarse_facetgid_entity=coarse_grid_labeling.d_to_dface_to_entity[Dc+1][tree]
          end 

          for poly_incident_face in poly_faces[poly_facet]
            if poly_incident_face == poly_facet
              ref_facetgid=cell_facets[ref_cell][gridap_facet]
              facet_to_entity[ref_facetgid]=coarse_facetgid_entity
              @debug "[FACE] info.tree_boundary=$(info.tree_boundary) facet_to_entity[$(ref_facetgid)]=$(coarse_facetgid_entity)"
            elseif (Dc==3 && poly_incident_face in Gridap.ReferenceFEs.get_dimrange(polytope,1))
              poly_first_edget=first(Gridap.ReferenceFEs.get_dimrange(polytope,1))
              edget=poly_incident_face-poly_first_edget+1
              ref_edgetgid=cell_edgets[ref_cell][edget]
              if (edget_to_entity[ref_edgetgid]==0)
                edget_to_entity[ref_edgetgid]=coarse_facetgid_entity
                @debug "[EDGE] info.tree_boundary=$(info.tree_boundary) edget_to_entity[$(ref_edgetgid)]=$(coarse_facetgid_entity)"
              end
            else
              ref_cornergid=cell_vertices[ref_cell][poly_incident_face]
              if (vertex_to_entity[ref_cornergid]==0)
                @debug "[CORNER ON FACE] info.tree_boundary=$(info.tree_boundary) vertex_to_entity[$(ref_cornergid)]=$(coarse_facetgid_entity)"
                vertex_to_entity[ref_cornergid]=coarse_facetgid_entity
              end
            end
          end
        end 

        nsides=info.sides.elem_count
        for iside=1:nsides
          if (iside==2 && 
            sides[iside].is_hanging == 0 &&
            sides[1].is_hanging == 0)
            break
          end
          face=sides[iside].face+1
          tree=sides[iside].treeid+1
          if (sides[iside].is_hanging == 0)
            data=sides[iside].is.full
            if data.is_ghost==1
              ref_cell=pXest.local_num_quadrants+data.quadid+1
            else
              ref_cell=owned_trees_offset[tree]+data.quadid+1
            end 
            @debug "nsides=$(nsides) sides[$(iside)].is_hanging == $(sides[iside].is_hanging) process_facet(tree=$(tree),face=$(face),ref_cell=$(ref_cell))"
            process_facet(tree,face,ref_cell,info)
          else 
            for i=1:length(sides[iside].is.hanging.quadid)
              quadid=sides[iside].is.hanging.quadid[i]
              if (sides[iside].is.hanging.is_ghost[i]==1)
                ref_cell=pXest.local_num_quadrants+quadid+1
              else 
                ref_cell=owned_trees_offset[tree]+quadid+1
              end
              @debug "nsides=$(nsides) sides[$(iside)].is_hanging == $(sides[iside].is_hanging) process_facet(tree=$(tree),face=$(face),ref_cell=$(ref_cell))"
              process_facet(tree,face,ref_cell,info)
            end
          end
        end 
        nothing
    end

    #  C-callable face callback
    cface_callback = @cfunction($jface_callback,
                                 Cvoid,
                                 (Ptr{p8est_iter_face_info_t},Ptr{Cvoid}))

    # Iterate over cells
    num_cells=Gridap.Geometry.num_faces(topology,Dc)
    cell_to_entity=zeros(Int,num_cells)

    # Cell iterator callback
    function jcell_callback(pinfo     :: Ptr{p8est_iter_volume_info_t},
                            user_data :: Ptr{Cvoid})
      info=pinfo[]
      tree=info.treeid+1
      cell=owned_trees_offset[tree]+info.quadid+1
      cell_to_entity[cell]=coarse_grid_labeling.d_to_dface_to_entity[Dc+1][tree]
      nothing
    end
    ccell_callback = @cfunction($jcell_callback,
                                 Cvoid,
                                 (Ptr{p8est_iter_volume_info_t},Ptr{Cvoid}))

    iterator_mode=Ref{Int}(ITERATOR_RESTRICT_TO_BOUNDARY)
    if (Dc==2)
       p4est_iterate(ptr_pXest,ptr_pXest_ghost,iterator_mode,C_NULL,cface_callback,C_NULL)
       p4est_iterate(ptr_pXest,ptr_pXest_ghost,C_NULL,ccell_callback,C_NULL,ccorner_callback)
       iterator_mode[]=ITERATOR_RESTRICT_TO_INTERIOR
       p4est_iterate(ptr_pXest,ptr_pXest_ghost,iterator_mode,C_NULL,cface_callback,C_NULL)
    else
       p8est_iterate(ptr_pXest,ptr_pXest_ghost,iterator_mode,C_NULL,cface_callback,C_NULL,C_NULL)
       p8est_iterate(ptr_pXest,ptr_pXest_ghost,iterator_mode,C_NULL,C_NULL,cedge_callback,C_NULL)
       p8est_iterate(ptr_pXest,ptr_pXest_ghost,C_NULL,ccell_callback,C_NULL,C_NULL,ccorner_callback)
       iterator_mode[]=ITERATOR_RESTRICT_TO_INTERIOR
       p8est_iterate(ptr_pXest,ptr_pXest_ghost,iterator_mode,C_NULL,cface_callback,C_NULL,C_NULL)
    end
    if (Dc==2)
      vertex_to_entity, facet_to_entity, cell_to_entity
    else
      vertex_to_entity, edget_to_entity, facet_to_entity, cell_to_entity
    end
  end

  vertex_to_entity  = map(x->x[1]   , faces_to_entity)
  if Dc==3
    edget_to_entity = map(x->x[2]   , faces_to_entity)
  end
  facet_to_entity   = map(x->x[Dc]  , faces_to_entity)
  cell_to_entity    = map(x->x[Dc+1], faces_to_entity)


  polytope = Dc==2 ? QUAD : HEX
  update_face_to_entity_with_ghost_data!(vertex_to_entity,
                                          cell_prange,
                                          num_faces(polytope,0),
                                          cell_to_faces(topology,Dc,0))
  if Dc==3
    update_face_to_entity_with_ghost_data!(edget_to_entity,
                                            cell_prange,
                                            num_faces(polytope,1),
                                            cell_to_faces(topology,Dc,1))   
  end

  update_face_to_entity_with_ghost_data!(facet_to_entity,
                                         cell_prange,
                                         num_faces(polytope,Dc-1),
                                         cell_to_faces(topology,Dc,Dc-1))
 
  update_face_to_entity_with_ghost_data!(cell_to_entity,
                                        cell_prange,
                                        num_faces(polytope,Dc),
                                        cell_to_faces(topology,Dc,Dc))

#  map(vertex_to_entity,facet_to_entity,cell_to_entity) do vertex_to_entity,facet_to_entity,cell_to_entity
#    @assert all(vertex_to_entity .!= 0)
#    @assert all(facet_to_entity  .!= 0)
#    @assert all(cell_to_entity   .!= 0)
#  end 

  if (Dc==2)
    faces_to_entity=[vertex_to_entity,facet_to_entity,cell_to_entity]
  else
    faces_to_entity=[vertex_to_entity,edget_to_entity,facet_to_entity,cell_to_entity]
  end

  face_labeling = map(faces_to_entity...) do faces_to_entity...
    d_to_dface_to_entity       = Vector{Vector{Int}}(undef,Dc+1)
    d_to_dface_to_entity[1]    = faces_to_entity[1]
    if (Dc==3)
      d_to_dface_to_entity[2]  = faces_to_entity[2]
    end
    d_to_dface_to_entity[Dc]   = faces_to_entity[Dc]
    d_to_dface_to_entity[Dc+1] = faces_to_entity[Dc+1]
    Gridap.Geometry.FaceLabeling(d_to_dface_to_entity,
                                 copy(coarse_grid_labeling.tag_to_entities),
                                 copy(coarse_grid_labeling.tag_to_name))
  end
  face_labeling
end

function cell_to_faces(topology,cell_dim,face_dim)
  map(topology) do topology
   Gridap.Geometry.get_faces(topology,cell_dim,face_dim)
  end
end

function _fill_data!(data,entry::Integer,k)
  data[k]=entry
  k=k+1
end

function _fill_data!(data,entry::Vector{<:Integer},k)
  for i=1:length(entry)
    data[k]=entry[i]
    k=k+1
  end
  k
end

function init_cell_to_face_entity(num_faces_x_cell,
                                  cell_to_faces,
                                  face_to_entity)
  ptrs = Vector{Int32}(undef,length(cell_to_faces) + 1)
  ptrs[2:end] .= num_faces_x_cell
  Gridap.Arrays.length_to_ptrs!(ptrs)
  data = Vector{eltype(face_to_entity)}(undef, ptrs[end] - 1)
  cell_to_face_entity = lazy_map(Broadcasting(Reindex(face_to_entity)),cell_to_faces)
  k = 1
  for i = 1:length(cell_to_face_entity)
    for j = 1:length(cell_to_face_entity[i])
      k=_fill_data!(data,cell_to_face_entity[i][j],k)
    end
  end
  return PArrays.JaggedArray(data, ptrs)
end

function update_face_to_entity!(face_to_entity, cell_to_faces, cell_to_face_entity)
  for cell in 1:length(cell_to_faces)
      i_to_entity = cell_to_face_entity[cell]
      pini = cell_to_faces.ptrs[cell]
      pend = cell_to_faces.ptrs[cell + 1] - 1
      for (i, p) in enumerate(pini:pend)
        lid = cell_to_faces.data[p]
        face_to_entity[lid] = i_to_entity[i]
      end
  end
end

function update_face_to_entity_with_ghost_data!(
   face_to_entity,cell_prange,num_faces_x_cell,cell_to_faces)

   part_to_cell_to_entity = map(init_cell_to_face_entity,
                                      map(x->num_faces_x_cell,partition(cell_prange)),
                                      cell_to_faces,
                                      face_to_entity)
   fetch_cache = fetch_vector_ghost_values_cache(part_to_cell_to_entity,
                                                 partition(cell_prange))
   fetch_vector_ghost_values!(part_to_cell_to_entity,fetch_cache) |> wait
   map(update_face_to_entity!,
             face_to_entity,
             cell_to_faces,
             part_to_cell_to_entity)
end

function setup_ptr_pXest_objects(pXest_type::PXestType,
                                 comm,
                                 coarse_discrete_model,
                                 num_uniform_refinements)
  ptr_pXest_connectivity=setup_pXest_connectivity(coarse_discrete_model)
  # Create a new forest
  ptr_pXest = setup_pXest(pXest_type,comm,ptr_pXest_connectivity,num_uniform_refinements)
  # Build the ghost layer
  ptr_pXest_ghost=setup_pXest_ghost(pXest_type,ptr_pXest)
  ptr_pXest_lnodes=setup_pXest_lnodes(pXest_type, ptr_pXest, ptr_pXest_ghost)
  ptr_pXest_connectivity, ptr_pXest, ptr_pXest_ghost, ptr_pXest_lnodes
end

function setup_distributed_discrete_model(pXest_type::PXestType,
                                          parts,
                                          coarse_discrete_model,
                                          ptr_pXest_connectivity,
                                          ptr_pXest,
                                          ptr_pXest_ghost,
                                          ptr_pXest_lnodes)
  
   cell_prange = setup_cell_prange(pXest_type,parts,ptr_pXest,ptr_pXest_ghost)
   cell_vertex_gids=generate_cell_vertex_gids(ptr_pXest_lnodes,cell_prange)
   cell_corner_lids = generate_cell_corner_lids(cell_vertex_gids)
   cell_vertex_coordinates = generate_cell_vertex_coordinates(pXest_type,
                                              cell_corner_lids,
                                              ptr_pXest_connectivity,
                                              ptr_pXest,
                                              ptr_pXest_ghost)

   grid,topology = generate_grid_and_topology(pXest_type,
                                              cell_corner_lids,
                                              cell_vertex_coordinates)

   face_labeling = generate_face_labeling(pXest_type,
                                          parts,
                                          cell_prange,
                                          coarse_discrete_model,
                                          topology,
                                          ptr_pXest,
                                          ptr_pXest_ghost)

   local_models = map(grid,topology,face_labeling) do grid, topology, face_labeling
     Gridap.Geometry.UnstructuredDiscreteModel(grid,topology,face_labeling)
   end
   return GridapDistributed.DistributedDiscreteModel(local_models,cell_prange)
end


"""
"""
function UniformlyRefinedForestOfOctreesDiscreteModel(
    parts::AbstractVector{<:Integer},
    coarse_discrete_model::DiscreteModel{Dc,Dp},
    num_uniform_refinements::Int
) where {Dc,Dp}

  comm = parts.comm

  pXest_type = _dim_to_pXest_type(Dc)

  ptr_pXest_connectivity,
     ptr_pXest,
       ptr_pXest_ghost,
         ptr_pXest_lnodes = setup_ptr_pXest_objects(pXest_type,
                                                    comm,
                                                    coarse_discrete_model,
                                                    num_uniform_refinements)

  # Write forest to VTK file
  # p4est_vtk_write_file(unitsquare_forest, C_NULL, "my_step")
  dmodel=setup_distributed_discrete_model(pXest_type,
                                          parts,
                                          coarse_discrete_model,
                                          ptr_pXest_connectivity,
                                          ptr_pXest,
                                          ptr_pXest_ghost,
                                          ptr_pXest_lnodes)

  pXest_lnodes_destroy(pXest_type,ptr_pXest_lnodes)
  pXest_ghost_destroy(pXest_type,ptr_pXest_ghost)
  pXest_destroy(pXest_type,ptr_pXest)
  pXest_connectivity_destroy(pXest_type,ptr_pXest_connectivity)
  dmodel
end
