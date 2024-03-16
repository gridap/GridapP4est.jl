abstract type PXestType end;
struct P4estType <: PXestType end ;
struct P6estType <: PXestType end;
struct P8estType <: PXestType end; 
const  P4P8estType = Union{P4estType,P8estType}     
const  P6P8estType = Union{P6estType,P8estType}

abstract type PXestRefinementRuleType end;
struct PXestUniformRefinementRuleType <: PXestRefinementRuleType end;
struct PXestVerticalRefinementRuleType <: PXestRefinementRuleType end;
struct PXestHorizontalRefinementRuleType <: PXestRefinementRuleType end;  

function pXest_destroy(pXest_type::P4estType, ptr_pXest)
  p4est_destroy(ptr_pXest)
end

function pXest_destroy(pXest_type::P6estType, ptr_pXest)
  p6est_destroy(ptr_pXest)
end

function pXest_destroy(pXest_type::P8estType, ptr_pXest)
  p8est_destroy(ptr_pXest)
end


function pXest_lnodes_destroy(pXest_type::P4estType, ptr_pXest_lnodes)
  p4est_lnodes_destroy(ptr_pXest_lnodes)
end

function pXest_lnodes_destroy(pXest_type::P6P8estType, ptr_pXest_lnodes)
  p8est_lnodes_destroy(ptr_pXest_lnodes)
end

function pXest_ghost_destroy(pXest_type::P4estType,ptr_pXest_ghost)
  p4est_ghost_destroy(ptr_pXest_ghost)
end 

function pXest_ghost_destroy(pXest_type::P6estType,ptr_pXest_ghost)
  p6est_ghost_destroy(ptr_pXest_ghost)
end

function pXest_ghost_destroy(pXest_type::P8estType,ptr_pXest_ghost)
  p8est_ghost_destroy(ptr_pXest_ghost)
end

function pXest_connectivity_destroy(pXest_type::P4estType, ptr_pXest_connectivity)
  p4est_connectivity_destroy(ptr_pXest_connectivity)
end

function pXest_connectivity_destroy(pXest_type::P6estType, ptr_pXest_connectivity)
  p6est_connectivity_destroy(ptr_pXest_connectivity)
end

function pXest_connectivity_destroy(pXest_type::P8estType, ptr_pXest_connectivity)
  p8est_connectivity_destroy(ptr_pXest_connectivity)
end

function setup_pXest(pXest_type::P4estType, comm, connectivity, num_uniform_refinements)
  p4est_new_ext(comm,
                connectivity,
                Cint(0), 
                Cint(num_uniform_refinements), 
                Cint(1), 
                Cint(0),
                C_NULL, 
                C_NULL)
end

function setup_pXest(pXest_type::P6estType, comm, connectivity, num_uniform_refinements)
    p6est_new_ext(comm,
                  connectivity,
                  Cint(0), 
                  Cint(num_uniform_refinements), # min_level 
                  Cint(num_uniform_refinements), # min_zlevel
                  Cint(1),                       # num_zroot
                  Cint(0),                       # fill_uniform
                  Cint(1),                       # data_size 
                  C_NULL,                        # init_fn
                  C_NULL)                        # user_pointer
end

function setup_pXest(pXest_type::P6estType, 
                    comm, connectivity, 
                    num_horizontal_uniform_refinements,
                    num_vertical_uniform_refinements)
    p6est_new_ext(comm,
                  connectivity,
                  Cint(0), 
                  Cint(num_horizontal_uniform_refinements), # min_level 
                  Cint(num_vertical_uniform_refinements),   # min_zlevel
                  Cint(1),                       # num_zroot
                  Cint(0),                       # fill_uniform
                  Cint(1),                       # data_size 
                  C_NULL,                        # init_fn
                  C_NULL)                        # user_pointer
end

function setup_pXest(pXest_type::P8estType, comm, connectivity, num_uniform_refinements)
    p8est_new_ext(comm,
                  connectivity,
                  Cint(0), Cint(num_uniform_refinements), Cint(1), Cint(0),
                  C_NULL, C_NULL)
end

function setup_pXest_ghost(pXest_type::P4estType, ptr_pXest)
  p4est_ghost_new(ptr_pXest,P4est_wrapper.P4EST_CONNECT_FULL)
end 

function setup_pXest_ghost(pXest_type::P6estType, ptr_pXest)
  p6est_ghost_new(ptr_pXest,P4est_wrapper.P4EST_CONNECT_FULL)
end

function setup_pXest_ghost(pXest_type::P8estType, ptr_pXest)
  p8est_ghost_new(ptr_pXest,P4est_wrapper.P8EST_CONNECT_FULL)
end 

function setup_pXest_lnodes(pXest_type::P4estType, ptr_pXest, ptr_pXest_ghost)
  p4est_lnodes_new(ptr_pXest, ptr_pXest_ghost, Cint(1))
end

function setup_pXest_lnodes(pXest_type::P8estType, ptr_pXest, ptr_pXest_ghost)
  p8est_lnodes_new(ptr_pXest, ptr_pXest_ghost, Cint(1))
end

function setup_pXest_lnodes_nonconforming(pXest_type::P4estType, ptr_pXest, ptr_pXest_ghost)
  p4est_lnodes_new(ptr_pXest, ptr_pXest_ghost, Cint(-2))
end 

function setup_pXest_lnodes_nonconforming(pXest_type::P6estType, ptr_pXest, ptr_pXest_ghost)
  p6est_lnodes_new(ptr_pXest, ptr_pXest_ghost, Cint(2))
end

function setup_pXest_lnodes_nonconforming(pXest_type::P8estType, ptr_pXest, ptr_pXest_ghost)
  p8est_lnodes_new(ptr_pXest, ptr_pXest_ghost, Cint(-3))
end 

function fill_tree_to_vertex!(conn,trian::Triangulation{D,D}) where D 
  cell_nodes_ids=Gridap.Geometry.get_cell_node_ids(trian)
  tree_to_vertex=unsafe_wrap(Array, conn.tree_to_vertex, length(cell_nodes_ids)*(2^D))
  c=Gridap.Arrays.array_cache(cell_nodes_ids)
  current=1
  for j=1:length(cell_nodes_ids)
    ids=Gridap.Arrays.getindex!(c,cell_nodes_ids,j)
    for id in ids
      tree_to_vertex[current]=p4est_topidx_t(id-1)
      current=current+1
    end
  end
end 

function fill_coordinates!(conn,trian::Triangulation{D,D}) where D 
  node_coordinates=Gridap.Geometry.get_node_coordinates(trian)
  vertices=unsafe_wrap(Array, conn.vertices, length(node_coordinates)*3)
  current=1
  for i=1:length(node_coordinates)
    p=node_coordinates[i]
    for j=1:D
      vertices[current]=Cdouble(p[j])
      current=current+1
    end
    if (D==2)
      vertices[current]=Cdouble(0.0) # Z coordinate always to 0.0 in 2D
      current=current+1
    end
  end
end 

function fill_tree_to_tree_and_to_face!(conn,trian::Triangulation{D,D}) where D
    # /*
    #  * Fill tree_to_tree and tree_to_face to make sure we have a valid
    #  * connectivity.
    #  */
    PXEST_FACES=2*D
    tree_to_tree=unsafe_wrap(Array, conn.tree_to_tree, conn.num_trees*PXEST_FACES )
    tree_to_face=unsafe_wrap(Array, conn.tree_to_face, conn.num_trees*PXEST_FACES )
    for tree=1:conn.num_trees
      for face=1:PXEST_FACES
        tree_to_tree[PXEST_FACES * (tree-1) + face] = tree-1
        tree_to_face[PXEST_FACES * (tree-1) + face] = face-1
      end
    end
end

function setup_pXest_connectivity(pXest_type::P4estType, coarse_discrete_model::DiscreteModel{2,2})
    trian=Triangulation(coarse_discrete_model)
    pconn=p4est_connectivity_new(
        p4est_topidx_t(num_nodes(trian)),                 # num_vertices
        p4est_topidx_t(num_cells(coarse_discrete_model)), # num_trees
        p4est_topidx_t(0),
        p4est_topidx_t(0))
    conn=pconn[]
    fill_tree_to_vertex!(conn, trian) 
    fill_coordinates!(conn, trian)
    fill_tree_to_tree_and_to_face!(conn, trian)
    p4est_connectivity_complete(pconn)
    @assert Bool(p4est_connectivity_is_valid(pconn))
    pconn 
end 

function setup_pXest_connectivity(pXest_type::P6estType, 
                                  coarse_discrete_model::DiscreteModel{2,2},
                                  extrusion_vector::Vector{Float64})
    @assert length(extrusion_vector)==3
    pconn4=setup_pXest_connectivity(P4estType(), coarse_discrete_model)
    p6est_connectivity_new(pconn4,C_NULL,extrusion_vector)
end 

function setup_pXest_connectivity(pXest_type::P8estType, coarse_discrete_model::DiscreteModel{3,3})
    trian=Triangulation(coarse_discrete_model)
    pconn=p8est_connectivity_new(
        p4est_topidx_t(length(node_coordinates)),         # num_vertices
        p4est_topidx_t(num_cells(coarse_discrete_model)), # num_trees
        p4est_topidx_t(0),
        p4est_topidx_t(0),
        p4est_topidx_t(0),
        p4est_topidx_t(0))
    conn=pconn[]
    fill_tree_to_vertex!(conn, trian) 
    fill_coordinates!(conn, trian)
    fill_tree_to_tree_and_to_face!(conn, trian)
    p8est_connectivity_complete(pconn)
    @assert Bool(p8est_connectivity_is_valid(pconn))
    pconn 
end

function pXest_reset_data!(::P4estType, ptr_pXest, data_size, init_fn_c, user_pointer)
  p4est_reset_data(ptr_pXest, data_size, init_fn_c, user_pointer)
end

function pXest_reset_data!(::P6estType, ptr_pXest, data_size, init_fn_c, user_pointer)
  p6est_reset_data(ptr_pXest, data_size, init_fn_c, user_pointer)
end

function pXest_reset_data!(::P8estType, ptr_pXest, data_size, init_fn_c, user_pointer)
  p8est_reset_data(ptr_pXest, data_size, init_fn_c, user_pointer)
end


function pXest_refine!(::P4estType, ptr_pXest, refine_fn_c, refine_replace_fn_c; init_fn_c=C_NULL) 
  p4est_refine_ext(ptr_pXest, Cint(0), Cint(-1), refine_fn_c, init_fn_c, refine_replace_fn_c)
end

function pXest_refine!(::P8estType, ptr_pXest, refine_fn_c, refine_replace_fn_c; init_fn_c=C_NULL) 
  p8est_refine_ext(ptr_pXest, Cint(0), Cint(-1), refine_fn_c, init_fn_c, refine_replace_fn_c)
end

function p6est_vertically_refine!(ptr_pXest, refine_fn_c, refine_replace_fn_c; init_fn_c=C_NULL) 
  p6est_refine_layers_ext(ptr_pXest, Cint(0), Cint(-1), refine_fn_c, init_fn_c, refine_replace_fn_c)
end

function p6est_horizontally_refine!(ptr_pXest, refine_fn_c, refine_replace_fn_c; init_fn_c=C_NULL) 
  p6est_refine_columns_ext(ptr_pXest, Cint(0), Cint(-1), refine_fn_c, init_fn_c, refine_replace_fn_c)
end

function p6est_vertically_coarsen!(ptr_pXest, coarsen_fn_c)
  p6est_coarsen_layers(ptr_pXest, Cint(0), coarsen_fn_c, C_NULL)
end

function p6est_horizontally_coarsen!(ptr_pXest, coarsen_fn_c)
  p6est_coarsen_columns(ptr_pXest, Cint(0), coarsen_fn_c, C_NULL)
end

function pXest_coarsen!(::P4estType, ptr_pXest, coarsen_fn_c) 
  p4est_coarsen(ptr_pXest, Cint(0), coarsen_fn_c, C_NULL)
end

function pXest_coarsen!(::P8estType, ptr_pXest, coarsen_fn_c)
  p8est_coarsen(ptr_pXest, Cint(0), coarsen_fn_c, C_NULL)
end 

function pXest_copy(::P4estType, ptr_pXest)
  p4est_copy(ptr_pXest, Cint(1))
end 

function pXest_copy(::P6estType, ptr_pXest)
  p6est_copy(ptr_pXest, Cint(1))
end

function pXest_copy(::P8estType, ptr_pXest)
  p8est_copy(ptr_pXest, Cint(1))
end

function pXest_balance!(::P4estType, ptr_pXest; k_2_1_balance=0)
  if (k_2_1_balance==0)
    p4est_balance(ptr_pXest, P4est_wrapper.P4EST_CONNECT_FULL, C_NULL) 
  else
    @assert k_2_1_balance==1 
    p4est_balance(ptr_pXest, P4est_wrapper.P4EST_CONNECT_FACE, C_NULL)
  end
end

function pXest_balance!(::P6estType, ptr_pXest; k_2_1_balance=0)
  @assert k_2_1_balance==0 or k_2_1_balance==2
  if (k_2_1_balance==0)
    p6est_balance(ptr_pXest, P4est_wrapper.P8EST_CONNECT_FULL, C_NULL) 
  elseif (k_2_1_balance==2)
    p6est_balance(ptr_pXest, P4est_wrapper.P8EST_CONNECT_FACE, C_NULL)
  end
end

function pXest_balance!(::P8estType, ptr_pXest; k_2_1_balance=0)
  if (k_2_1_balance==0)
    p8est_balance(ptr_pXest, P4est_wrapper.P8EST_CONNECT_FULL, C_NULL) 
  elseif (k_2_1_balance==1)
    p8est_balance(ptr_pXest, P4est_wrapper.P8EST_CONNECT_EDGE, C_NULL)
  else 
    @assert k_2_1_balance==2
    p8est_balance(ptr_pXest, P4est_wrapper.P8EST_CONNECT_FACE, C_NULL)  
  end
end

function pXest_partition!(::P4estType, ptr_pXest)
  p4est_partition(ptr_pXest, 0, C_NULL)
end 

function pXest_partition!(::P6estType, ptr_pXest)
  p6est_partition(ptr_pXest, C_NULL)
end 

function pXest_partition!(::P8estType, ptr_pXest)
  p8est_partition(ptr_pXest, 0, C_NULL)
end 


function pXest_reset_callbacks(::P4estType)
  # Variables which are updated accross calls to init_fn_callback
  current_quadrant_index_within_tree = Cint(0)
  current_quadrant_index_among_trees = Cint(0)

  # This C callback function is called once per quadtree quadrant. Here we are assuming
  # that p4est->user_pointer has been set prior to the first call to this call
  # back function to an array of ints with as many entries as forest quadrants. This call back function
  # initializes the quadrant->p.user_data void * pointer of all quadrants such that it
  # points to the corresponding entry in the global array mentioned in the previous sentence
  function init_fn_callback(forest_ptr::Ptr{p4est_t},
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
  @cfunction($init_fn_callback, 
             Cvoid, (Ptr{p4est_t}, p4est_topidx_t, Ptr{p4est_quadrant_t}))
end

function pXest_reset_callbacks(::P8estType)
  # Variables which are updated accross calls to init_fn_callback
  current_quadrant_index_within_tree = Cint(0)
  current_quadrant_index_among_trees = Cint(0)

  function init_fn_callback(forest_ptr::Ptr{p8est_t},
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
  init_fn_callback_c = @cfunction($init_fn_callback, 
                                  Cvoid, (Ptr{p8est_t}, p4est_topidx_t, Ptr{p8est_quadrant_t}))
end 


function p2est_quadrant_is_equal(a,b)
  a[].z==b[].z && a[].level==b[].level
end 

function P4EST_QUADRANT_LEN(l)
  p4est_qcoord_t(1) << (P4est_wrapper.P4EST_MAXLEVEL-l)
end 

function p2est_quadrant_is_ancestor(a,b)
  if (a[].level>=b[].level)
    return false
  end 
  return (b[].z >= a[].z && 
          b[].z < a[].z + P4EST_QUADRANT_LEN(a[].level))
end

function p6est_vertically_adapt_reset_callbacks()
  current_quadrant_index_among_trees  = Cint(-1)
  current_quadrant_index_within_tree  = Cint(0)
  current_layer_within_column         = Cint(0)
  current_layer                       = Cint(0)
  previous_quadrant                   = Ref{p4est_quadrant_t}()
  

  function init_fn_callback(forest_ptr::Ptr{p6est_t},
    which_tree::p4est_topidx_t,
    column_ptr::Ptr{p4est_quadrant_t},
    layer_ptr::Ptr{p2est_quadrant_t})

    forest = forest_ptr[]
    columns = forest.columns[]
    tree = p4est_tree_array_index(columns.trees, which_tree)[]
    quadrant = column_ptr[]
    layer = layer_ptr[]
  
    if (current_quadrant_index_among_trees==-1 || 
          p4est_quadrant_compare(previous_quadrant,column_ptr) != 0)
  
      previous_quadrant = column_ptr
      current_quadrant_index_among_trees = current_quadrant_index_among_trees+1
      current_quadrant_index_within_tree = (current_quadrant_index_within_tree + 1) % (tree.quadrants.elem_count)
      q = p4est_quadrant_array_index(tree.quadrants, current_quadrant_index_within_tree)
      #@assert p4est_quadrant_compare(q,column_ptr) == 0
      current_layer_within_column = 0
    end
  
    q = p4est_quadrant_array_index(tree.quadrants, current_quadrant_index_within_tree)
    #@assert p4est_quadrant_compare(q,column_ptr) == 0


    user_data  = unsafe_wrap(Array, Ptr{Cint}(forest.user_pointer), current_layer+1)[current_layer+1]

  
    f,l=P6EST_COLUMN_GET_RANGE(column_ptr[])
    q2_ptr=p2est_quadrant_array_index(forest.layers[], f+current_layer_within_column)
    @assert p2est_quadrant_is_equal(q2_ptr,layer_ptr)

    unsafe_store!(Ptr{Cint}(q2_ptr[].p.user_data), user_data, 1)
  
  
    current_layer_within_column=current_layer_within_column+1
    current_layer=current_layer+1                          
        
    if ((which_tree+1)==columns.connectivity[].num_trees && 
        current_quadrant_index_within_tree==tree.quadrants.elem_count && 
        current_layer_within_column==l-f)
      current_quadrant_index_among_trees = Cint(-1)
      current_layer=Cint(0)
    end 
  
    return nothing
  end

  @cfunction($init_fn_callback, Cvoid, 
             (Ptr{p6est_t}, p4est_topidx_t, Ptr{p4est_quadrant_t},Ptr{p2est_quadrant_t}))

end 


function p6est_horizontally_adapt_reset_callbacks()
  current_quadrant_index_among_trees  = Cint(-1)
  current_quadrant_index_within_tree  = Cint(0)
  current_layer_within_column         = Cint(0)
  current_layer                       = Cint(0)
  previous_quadrant                   = Ref{p4est_quadrant_t}()
  

  function init_fn_callback(forest_ptr::Ptr{p6est_t},
    which_tree::p4est_topidx_t,
    column_ptr::Ptr{p4est_quadrant_t},
    layer_ptr::Ptr{p2est_quadrant_t})

    forest = forest_ptr[]
    columns = forest.columns[]
    tree = p4est_tree_array_index(columns.trees, which_tree)[]
    quadrant = column_ptr[]
    layer = layer_ptr[]
  
    if (current_quadrant_index_among_trees==-1 || 
          p4est_quadrant_compare(previous_quadrant,column_ptr) != 0)
  
      previous_quadrant = column_ptr
      current_quadrant_index_among_trees = current_quadrant_index_among_trees+1
      current_quadrant_index_within_tree = (current_quadrant_index_within_tree + 1) % (tree.quadrants.elem_count)
      q = p4est_quadrant_array_index(tree.quadrants, current_quadrant_index_within_tree)
      current_layer_within_column = 0
    end
    
    f,l=P6EST_COLUMN_GET_RANGE(column_ptr[])
    if (current_layer_within_column==0)
       user_data = 
          unsafe_wrap(Array, Ptr{Cint}(forest.user_pointer), 
             current_quadrant_index_among_trees+1)[current_quadrant_index_among_trees+1]
       # We will use the first layer of the first column 
       # to decide whether we refine the column or not   
       q2_ptr=p2est_quadrant_array_index(forest.layers[], f)
       @assert p2est_quadrant_is_equal(q2_ptr,layer_ptr)
        
       unsafe_store!(Ptr{Cint}(q2_ptr[].p.user_data), user_data, 1)
    end
  
    current_layer_within_column=current_layer_within_column+1
    current_layer=current_layer+1                          
        
    if ((which_tree+1)==columns.connectivity[].num_trees && 
        current_quadrant_index_within_tree==tree.quadrants.elem_count && 
        current_layer_within_column==l-f)
      current_quadrant_index_among_trees = Cint(-1)
      current_layer=Cint(0)
      current_quadrant_index_within_tree = Cint(0)
    end 
    return nothing
  end
  @cfunction($init_fn_callback, Cvoid, 
             (Ptr{p6est_t}, p4est_topidx_t, Ptr{p4est_quadrant_t},Ptr{p2est_quadrant_t}))
end 


function p6est_vertically_refine_callbacks()
  function refine_layer_callback(p6est::Ptr{p6est_t},
    which_tree::p4est_topidx_t,
    column::Ptr{p4est_quadrant_t},
    layer::Ptr{p2est_quadrant_t})
    Cint(unsafe_wrap(Array, Ptr{Cint}(layer[].p.user_data), 1)[] == refine_flag)
  end
  
  refine_layer_fn_callback_c = @cfunction($refine_layer_callback, Cint, 
                                    (Ptr{p6est_t}, p4est_topidx_t, Ptr{p4est_quadrant_t}, Ptr{p2est_quadrant_t}))

  function refine_layer_replace_callback(p6est::Ptr{p6est_t},
    which_tree::p4est_topidx_t,
    num_outcolumns::Cint,
    num_outlayers::Cint,
    outcolumns::Ptr{Ptr{p4est_quadrant_t}},
    outlayers::Ptr{Ptr{p2est_quadrant_t}},
    num_incolumns::Cint,
    num_inlayers::Cint,
    incolumns::Ptr{Ptr{p4est_quadrant_t}},
    inlayers::Ptr{Ptr{p2est_quadrant_t}})

    @assert num_outcolumns==1
    @assert num_outlayers==1
    @assert num_incolumns==1
    @assert num_inlayers==2

    inlayers_array=unsafe_wrap(Array, inlayers, num_inlayers)
    for i=1:num_inlayers
      quadrant = inlayers_array[i][]
      if (quadrant.p.user_data) != C_NULL
          unsafe_store!(Ptr{Cint}(quadrant.p.user_data), nothing_flag, 1)
      end
    end
  end

  refine_layer_replace_callback_c = @cfunction($refine_layer_replace_callback, Cvoid, 
      (Ptr{p6est_t}, p4est_topidx_t, Cint, Cint, 
       Ptr{Ptr{p4est_quadrant_t}}, Ptr{Ptr{p2est_quadrant_t}}, Cint, Cint, 
       Ptr{Ptr{p4est_quadrant_t}}, Ptr{Ptr{p2est_quadrant_t}}))

  refine_layer_fn_callback_c, refine_layer_replace_callback_c
end

function p6est_horizontally_refine_callbacks()
  function refine_column_callback(p6est::Ptr{p6est_t},
    which_tree::p4est_topidx_t,
    column::Ptr{p4est_quadrant_t})
    forest=p6est[]
    f,l=P6EST_COLUMN_GET_RANGE(column[])
    q2_ptr=p2est_quadrant_array_index(forest.layers[], f)
    Cint(unsafe_wrap(Array, Ptr{Cint}(q2_ptr[].p.user_data), 1)[] == refine_flag)
  end

  refine_column_fn_callback_c = @cfunction($refine_column_callback, Cint, 
                                           (Ptr{p6est_t}, p4est_topidx_t, Ptr{p4est_quadrant_t}))

  function refine_column_replace_callback(p6est::Ptr{p6est_t},
    which_tree::p4est_topidx_t,
    num_outcolumns::Cint,
    num_outlayers::Cint,
    outcolumns::Ptr{Ptr{p4est_quadrant_t}},
    outlayers::Ptr{Ptr{p2est_quadrant_t}},
    num_incolumns::Cint,
    num_inlayers::Cint,
    incolumns::Ptr{Ptr{p4est_quadrant_t}},
    inlayers::Ptr{Ptr{p2est_quadrant_t}})

    @assert num_outcolumns==1
    @assert num_incolumns==4

    inlayers_array=unsafe_wrap(Array, inlayers, num_inlayers)
    for i=1:num_inlayers
      quadrant = inlayers_array[i][]
      if (quadrant.p.user_data) != C_NULL
          unsafe_store!(Ptr{Cint}(quadrant.p.user_data), nothing_flag, 1)
      end
    end
  end

  refine_column_replace_callback_c = @cfunction($refine_column_replace_callback, Cvoid, 
      (Ptr{p6est_t}, p4est_topidx_t, Cint, Cint, 
       Ptr{Ptr{p4est_quadrant_t}}, Ptr{Ptr{p2est_quadrant_t}}, Cint, Cint, 
       Ptr{Ptr{p4est_quadrant_t}}, Ptr{Ptr{p2est_quadrant_t}}))

  refine_column_fn_callback_c, refine_column_replace_callback_c
end 


function p6est_vertically_coarsen_callbacks()
  function coarsen_layer_callback(p6est::Ptr{p6est_t},
    which_tree::p4est_topidx_t,
    column::Ptr{p4est_quadrant_t},
    layer::Ptr{Ptr{p2est_quadrant_t}})

    num_children=2
    layers=unsafe_wrap(Array, layer, num_children)

    coarsen=Cint(1)
    for quadrant_index=1:num_children
       quadrant = layers[quadrant_index][]
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
  coarsen_fn_callback_c = @cfunction($coarsen_layer_callback, 
                                     Cint, 
                                     (Ptr{p6est_t}, p4est_topidx_t, Ptr{p4est_quadrant_t}, Ptr{Ptr{p2est_quadrant_t}}))

  return coarsen_fn_callback_c
end 


function p6est_horizontally_coarsen_callbacks()
  function coarsen_column_callback(p6est::Ptr{p6est_t},
    which_tree::p4est_topidx_t,
    columns::Ptr{Ptr{p4est_quadrant_t}})

    forest=p6est[]
    num_children=4
    columns_array=unsafe_wrap(Array, columns, num_children)

    coarsen=Cint(1)
    for quadrant_index=1:num_children
       column = columns_array[quadrant_index][]
       f,l = P6EST_COLUMN_GET_RANGE(column)
       q2_ptr = p2est_quadrant_array_index(forest.layers[], f)
       if (q2_ptr[].p.user_data) == C_NULL
        return Cint(0)
       end
       is_coarsen_flag = (unsafe_wrap(Array,Ptr{Cint}(q2_ptr[].p.user_data),1)[])==coarsen_flag
       if (!is_coarsen_flag) 
         return Cint(0)
       end
    end
    return coarsen
  end
  coarsen_fn_callback_c = @cfunction($coarsen_column_callback, 
                                     Cint, 
                                     (Ptr{p6est_t}, p4est_topidx_t, Ptr{Ptr{p4est_quadrant_t}}))

  return coarsen_fn_callback_c
end 


function pXest_coarsen_callbacks(::P4estType)
  function coarsen_callback(forest_ptr::Ptr{p4est_t},
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
  coarsen_fn_callback_c = @cfunction($coarsen_callback, 
                                     Cint, (Ptr{p4est_t}, p4est_topidx_t, Ptr{Ptr{p4est_quadrant_t}}))
end 

function pXest_coarsen_callbacks(::P8estType)
  function coarsen_callback(forest_ptr::Ptr{p8est_t},
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
  coarsen_fn_callback_c = @cfunction($coarsen_callback, 
                                    Cint, (Ptr{p8est_t}, p4est_topidx_t, Ptr{Ptr{p8est_quadrant_t}}))
end 


function pXest_refine_callbacks(::P4estType)
  function refine_callback(::Ptr{p4est_t},
    which_tree::p4est_topidx_t,
    quadrant_ptr::Ptr{p4est_quadrant_t})
    quadrant = quadrant_ptr[]
    return Cint(unsafe_wrap(Array, Ptr{Cint}(quadrant.p.user_data), 1)[] == refine_flag)
  end
  refine_callback_c = @cfunction($refine_callback, Cint, (Ptr{p4est_t}, p4est_topidx_t, Ptr{p4est_quadrant_t}))

  function refine_replace_callback(::Ptr{p4est_t},
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
  refine_replace_callback_c = 
    @cfunction($refine_replace_callback, Cvoid, (Ptr{p4est_t}, 
                                                    p4est_topidx_t, 
                                                    Cint, 
                                                    Ptr{Ptr{p4est_quadrant_t}}, 
                                                    Cint, 
                                                    Ptr{Ptr{p4est_quadrant_t}}))
  refine_callback_c, refine_replace_callback_c
end

function pXest_refine_callbacks(::P8estType)
  function refine_callback(::Ptr{p8est_t},
    which_tree::p4est_topidx_t,
    quadrant_ptr::Ptr{p8est_quadrant_t})
    quadrant = quadrant_ptr[]
    return Cint(unsafe_wrap(Array, Ptr{Cint}(quadrant.p.user_data), 1)[] == refine_flag)
  end
  refine_callback_c = @cfunction($refine_callback, Cint, (Ptr{p8est_t}, p4est_topidx_t, Ptr{p8est_quadrant_t}))

  function refine_replace_callback(::Ptr{p8est_t},
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

  refine_replace_callback_c = 
    @cfunction($refine_replace_callback, Cvoid, (Ptr{p8est_t}, 
                                                  p4est_topidx_t, 
                                                  Cint, 
                                                  Ptr{Ptr{p8est_quadrant_t}}, 
                                                  Cint, 
                                                  Ptr{Ptr{p8est_quadrant_t}}))

  refine_callback_c, refine_replace_callback_c
end 

function _unwrap_ghost_quadrants(::P4estType, pXest_ghost)
  Ptr{p4est_quadrant_t}(pXest_ghost.ghosts.array)
end

function _unwrap_ghost_quadrants(::P6estType, pXest_ghost)
  Ptr{p2est_quadrant_t}(pXest_ghost.ghosts.array)
end

function _unwrap_ghost_quadrants(::P8estType, pXest_ghost)
  Ptr{p8est_quadrant_t}(pXest_ghost.ghosts.array)
end

function _unwrap_global_first_quadrant(::P4P8estType, pXest)
  unsafe_wrap(Array,
              pXest.global_first_quadrant,
              pXest.mpisize+1)
end

function _unwrap_global_first_quadrant(::P6estType, pXest)
  unsafe_wrap(Array,
              pXest.global_first_layer,
              pXest.mpisize+1)
end


function setup_cell_prange(pXest_type::PXestType,
                           parts::AbstractVector{<:Integer},
                           ptr_pXest,
                           ptr_pXest_ghost)
  comm = parts.comm

  pXest_ghost = ptr_pXest_ghost[]
  pXest       = ptr_pXest[]

  # Obtain ghost quadrants
  ptr_ghost_quadrants = _unwrap_ghost_quadrants(pXest_type,pXest_ghost)
  proc_offsets = unsafe_wrap(Array, pXest_ghost.proc_offsets, pXest_ghost.mpisize+1)

  global_first_quadrant = _unwrap_global_first_quadrant(pXest_type,pXest)

  noids,firstgid,gho_to_glo,gho_to_own=map(parts) do part
    gho_to_glo = Vector{Int}(undef, pXest_ghost.ghosts.elem_count)
    gho_to_own = Vector{Int32}(undef, pXest_ghost.ghosts.elem_count)
    k=1
    for i=1:pXest_ghost.mpisize
      for j=proc_offsets[i]:proc_offsets[i+1]-1
        quadrant       = ptr_ghost_quadrants[j+1]
        piggy3         = quadrant.p.piggy3
        gho_to_glo[k]  = global_first_quadrant[i]+piggy3.local_num+1
        gho_to_own[k] = Int32(i)
        k=k+1
      end
    end
    global_first_quadrant[part+1]-global_first_quadrant[part],global_first_quadrant[part]+1,gho_to_glo,gho_to_own
  end |> tuple_of_arrays
  ngids = global_first_quadrant[end]

  partition=map(parts,noids,firstgid,gho_to_glo,gho_to_own) do part, noids, firstgid, gho_to_glo, gho_to_own
    owner = part
    own_indices=OwnIndices(ngids,owner,(collect(firstgid:firstgid+noids-1)))
    ghost_indices=GhostIndices(ngids,gho_to_glo,gho_to_own)
    OwnAndGhostIndices(own_indices,ghost_indices)
  end
  # This is required to provide the hint that the communication 
  # pattern underlying partition is symmetric, so that we do not have 
  # to execute the algorithm the reconstructs the reciprocal in the 
  # communication graph
  assembly_neighbors(partition;symmetric=true)
  PRange(partition)
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

function p6est_lnodes_decode(face_code,
                             hanging_face,
                             hanging_edge)
    @assert face_code >= 0
    if (face_code != 0)
        fc4  = face_code & 0x000f
        h    = Int16((face_code & 0x0010) >> 4)
        work = Int16(face_code >> 5)
        hanging_face .= -1
        hanging_edge .= -1
        p4est_lnodes_decode(fc4, view(hanging_face,3:length(hanging_face)))
        for f=0:3
            hf = hanging_face[f + 3]
            w = work & 0x0001
            if (hf >= 0)
              hanging_edge[p8est_face_edges[f+3,3]+1] = 2 + hf
              hanging_edge[p8est_face_edges[f+3,4]+1] = 2 + hf
              hanging_edge[p8est_face_edges[f+3,1⊻hf+1]+1] = 4
              if (w!=0)
                hanging_edge[p8est_face_edges[f+3,3⊻h+1]+1] = 4
                hanging_edge[p8est_face_edges[f+3,1⊻hf+1]+1] = 4
                hanging_edge[p8est_face_edges[f+3,hf+1]+1] = 2+h
                hanging_face[f + 3] = (hf << 1) | h
              else
                hanging_face[f + 3] = 4 + hf
              end
            elseif (w!=0)
              hanging_edge[p8est_face_edges[f+3,3⊻h+1]+1] = 4;
              hanging_edge[p8est_face_edges[f+3,1]+1] =
                max(hanging_edge[p8est_face_edges[f+3,1]+1], 2+h);
              hanging_edge[p8est_face_edges[f+3,2]+1] =
                max(hanging_edge[p8est_face_edges[f+3,2]+1], 2+h);
              hanging_face[f + 3] = 6 + h;
            end       
            work >>= 1;
        end

        for e=0:3
           if ((work & 0x0001)!=0)
             if (hanging_edge[e+1] < 0)
                hanging_edge[e+1] = h
             else
                @assert (hanging_edge[e+1] == 2 + h || hanging_edge[e+1] == 4)
             end
           end
           work >>= 1; 
        end 
        return 1
    else
        return 0
    end
end 

function pXest_lnodes_decode(::P6estType,face_code, hanging_face, hanging_edge)
  p6est_lnodes_decode(face_code, hanging_face, hanging_edge)
end

function pXest_lnodes_decode(::P8estType,face_code, hanging_face, hanging_edge)
  p8est_lnodes_decode(face_code, hanging_face, hanging_edge)
end

function num_cell_dims(::P4estType)
  2
end 

function num_cell_dims(::P6P8estType)
  3
end 

function _faces_to_cell_element_nodes(pXest_type::P4P8estType)
  Dc = num_cell_dims(pXest_type)
  nf = num_cell_faces(Val{Dc})
  collect(1:nf)
end 

function _faces_to_cell_element_nodes(pXest_type::P6estType)
  [13,15, 11,17, 5,23]
end 

function _edges_to_cell_element_nodes(pXest_type::P8estType)
  Dc = num_cell_dims(pXest_type)
  nf = num_cell_faces(Val{Dc})
  ne = num_cell_edges(Val{Dc})
  collect(nf+1:nf+ne)
end

function _edges_to_cell_element_nodes(pXest_type::P6estType)
  [2,8,20,26, 4,6,22,24, 10,12,16,18]
end

function _vertices_to_cell_element_nodes(pXest_type::P4P8estType)
  Dc = num_cell_dims(pXest_type)
  nf = num_cell_faces(Val{Dc})
  ne = num_cell_edges(Val{Dc})
  nv = num_cell_vertices(Val{Dc})
  collect(nf+ne+1:nf+ne+nv)
end

function _vertices_to_cell_element_nodes(pXest_type::P6estType)
  [1,3,7,9, 19,21,25,27]
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


function _build_map_from_faces_edges_to_cell_lface_ledge(pXest_type, vnodes, element_nodes, face_code)
  
  faces_to_cell_element_nodes = _faces_to_cell_element_nodes(pXest_type)
  edges_to_cell_element_nodes = _edges_to_cell_element_nodes(pXest_type)

  
  n_cell_faces    = num_cell_faces(Val{3})
  n_cell_edges    = num_cell_edges(Val{3})

  hanging_face = Vector{Cint}(undef, n_cell_faces)
  hanging_edge = Vector{Cint}(undef, n_cell_edges)

  # Build a map from faces to (cell,lface)
  p4est_gface_to_gcell_p4est_lface = Dict{Int,Tuple{Int,Int}}()
  p4est_gedge_to_gcell_p4est_ledge = Dict{Int,Tuple{Int,Int}}()
  for cell = 1:length(face_code)
    s = (cell - 1) * vnodes + 1
    e = cell * vnodes
    p4est_cell_nodes = view(element_nodes, s:e)
    p4est_cell_faces = view(p4est_cell_nodes,faces_to_cell_element_nodes)
    p4est_cell_edges = view(p4est_cell_nodes,edges_to_cell_element_nodes)
  
    
    has_hanging = pXest_lnodes_decode(pXest_type, face_code[cell], hanging_face, hanging_edge)
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


function subface_to_hanging_edges_within_subface(::P6estType)
  p6est_subface_to_hanging_edges_within_subface
end 

function subface_to_hanging_edges_within_subface(::P8estType)
  p8est_subface_to_hanging_edges_within_subface
end 

function subface_to_hanging_edges_within_face(::P6estType)
  p6est_subface_to_hanging_edges_within_face
end 

function subface_to_hanging_edges_within_face(::P8estType)
  p8est_subface_to_hanging_edges_within_face
end

const p6est_half_to_regular_vertices = [ 0 1; 2 3; 0 2; 1 3]

function generate_cell_faces_and_non_conforming_glue(pXest_type::PXestType, 
                                                     ptr_pXest_lnodes, 
                                                     cell_prange) 
  
  Dc = num_cell_dims(pXest_type)
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
    @debug "ENDES_WO_GHOSTS[$(part_id(indices))]: $(element_nodes)"
    @debug "ENDES_WITH_GHOSTS[$(part_id(indices))]: $(element_nodes_with_ghosts.data)"
    @debug "FCODS_WO_GHOSTS[$(part_id(indices))]: $(face_code)"
    @debug "FCODS_WITH_GHOSTS[$(part_id(indices))]: $(face_code_with_ghosts)"
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
           _build_map_from_faces_edges_to_cell_lface_ledge(pXest_type,
                                                           lnodes.vnodes, 
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
    
    faces_to_cell_element_nodes    = _faces_to_cell_element_nodes(pXest_type)
    if (Dc==3)
      edges_to_cell_element_nodes    = _edges_to_cell_element_nodes(pXest_type)
    end
    
    vertices_to_cell_element_nodes = _vertices_to_cell_element_nodes(pXest_type)

    for cell = 1:num_local_elements
      start_gridap_vertices = (cell - 1) * n_cell_vertices
      start_gridap_faces    = (cell - 1) * n_cell_faces

      s = (cell-1)*lnodes.vnodes + 1
      e = cell*lnodes.vnodes
      p4est_cell_nodes = view(element_nodes_with_ghosts.data, s:e)

      p4est_cell_faces = view(p4est_cell_nodes,faces_to_cell_element_nodes)
      p4est_cell_vertices = view(p4est_cell_nodes,vertices_to_cell_element_nodes)
    
      gridap_cell_vertices = view(gridap_cell_vertices_data,
        start_gridap_vertices+1:start_gridap_vertices+n_cell_vertices)
      gridap_cell_faces = view(gridap_cell_faces_data,
        start_gridap_faces+1:start_gridap_faces+n_cell_faces)

      if (Dc==2)  
        has_hanging = p4est_lnodes_decode(face_code_with_ghosts[cell], hanging_face)
      else
        has_hanging = pXest_lnodes_decode(pXest_type, face_code_with_ghosts[cell], hanging_face, hanging_edge)
        start_gridap_edges = (cell-1)*n_cell_edges
        gridap_cell_edges  = view(gridap_cell_edges_data, start_gridap_edges+1:start_gridap_edges+n_cell_edges)
        p4est_cell_edges   = view(p4est_cell_nodes,edges_to_cell_element_nodes)
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
        if (isa(pXest_type,P4P8estType))
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
        end

        if (Dc==3)
          for (p4est_ledge, half) in enumerate(hanging_edge)
            if (half != -1 && half !=4)
              hanging_vertex_lvertex_within_edge = hanging_lvertex_within_edge(half)
              p4est_lvertex = p8est_edge_corners[p4est_ledge,
                                                 hanging_vertex_lvertex_within_edge+1]
              gridap_cell_vertices[PXEST_2_GRIDAP_VERTEX[p4est_lvertex+1]] = hanging_vertex_code

              @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))] cell=$(cell) hanging p4est_ledge=$(p4est_ledge) half=$(half) hanging p4est_lvertex=$(PXEST_2_GRIDAP_VERTEX[p4est_lvertex+1])"
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
            owner_face = p4est_cell_faces[p4est_lface]

            if half in 0:3
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
              hanging_vertices_pairs_to_owner_face[(cell, PXEST_2_GRIDAP_VERTEX[p4est_hanging_lvertex+1])] = owner_face

              # Process hanging face
              hanging_faces_pairs_to_owner_face[(cell, PXEST_2_GRIDAP_FACE[p4est_lface])] = (owner_face,half+1)
            else
              # Anisotropic refinement
              @assert half in 4:7
              for regular_vertex_lvertex_within_face in p6est_half_to_regular_vertices[mod(half,4)+1,:]
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
              end
              subface = mod(half,2)+1

              # Process hanging face
              hanging_faces_pairs_to_owner_face[(cell, PXEST_2_GRIDAP_FACE[p4est_lface])] = (owner_face,subface)
            end  
           
            if (Dc==3)
              _subface_to_hanging_edges_within_subface = subface_to_hanging_edges_within_subface(pXest_type)
              _subface_to_hanging_edges_within_face    = subface_to_hanging_edges_within_face(pXest_type)

              for (i,ledge_within_face) in enumerate(_subface_to_hanging_edges_within_subface[mod(half,4)+1,:])
                p4est_ledge=p8est_face_edges[p4est_lface,ledge_within_face+1]

                @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))]  cell=$(cell) p4est_lface=$(p4est_lface) half=$(half) index=$(mod(half,4)+1) p4est_ledge=$(p4est_ledge) ledge_within_face=$(ledge_within_face) "

                gridap_ledge = PXEST_2_GRIDAP_EDGE[p4est_ledge+1]
                # Identify the two edges which are hanging within the face
                hanging_edges_cell_ledge_to_owner_face_half[(cell, gridap_ledge)] =
                    (owner_face,-_subface_to_hanging_edges_within_face[mod(half,4)+1,i])
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
              @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))]  cell=$(cell) hanging_edge=$(hanging_edge) p4est_ledge=$(p4est_ledge) gridap_cell_edges[PXEST_2_GRIDAP_EDGE[p4est_ledge]]: $(gridap_cell_edges[PXEST_2_GRIDAP_EDGE[p4est_ledge]]) half=$(half)"

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

                @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))]  cell=$(cell) hanging_edge=$(hanging_edge) p4est_ledge=$(p4est_ledge) gridap_cell_edges[PXEST_2_GRIDAP_EDGE[p4est_ledge]]: $(gridap_cell_edges[PXEST_2_GRIDAP_EDGE[p4est_ledge]]) half=$(half) owner_edge_subedge_pair=$(owner_edge_subedge_pair)"

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
      # The following loop needs gridap cell vertices to be already completed
      for key in keys(hanging_edges_cell_ledge_to_owner_face_half)
          (cell, ledge) = key
          (owner_p4est_gface_or_gedge, half) = hanging_edges_cell_ledge_to_owner_face_half[key]
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
              end
            end
          end
      end
      @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))] gridap_cell_edges_data: $(gridap_cell_edges_data)"
      
      for key in keys(hanging_edges_cell_ledge_to_owner_face_half)
        (cell, ledge) = key
        (owner_p4est_gface_or_gedge, half) = hanging_edges_cell_ledge_to_owner_face_half[key]
        @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))] own_length=$(own_length(indices)) cell=$(cell) ledge=$(ledge) owner_p4est_gface_or_gedge=$(owner_p4est_gface_or_gedge) half=$(half)"
        if (half>0) # hanging edge is within a coarser face 
          @assert half==1 || half==2
          owner_p4est_gedge = owner_p4est_gface_or_gedge
          owner_gedge_pair = (owner_p4est_gedge,half)
          if (haskey(owner_edge_subedge_to_cell_ledge,owner_gedge_pair))
            (owner_cell, owner_cell_ledge) = owner_edge_subedge_to_cell_ledge[owner_gedge_pair]
            if (owner_cell!=cell)
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

 function pXest_tree_array_index(::P4estType, pXest, i)
  p4est_tree_array_index(pXest.trees, i)
end 

function pXest_tree_array_index(::P8estType, pXest, i)
  p8est_tree_array_index(pXest.trees, i)
end

function pXest_tree_array_index(::P6estType, pXest, i)
  p4est_tree_array_index(pXest.columns[].trees, i)
end

function pXest_quadrant_array_index(::P4estType, tree, i)
  p4est_quadrant_array_index(tree.quadrants, i)
end 

function pXest_quadrant_array_index(::P6estType, tree, i)
  p4est_quadrant_array_index(tree.quadrants, i)
end

function pXest_quadrant_array_index(::P8estType, tree, i)
  p8est_quadrant_array_index(tree.quadrants, i)
end

function pXest_quadrant_is_parent(::P4estType, q1, q2)
  p4est_quadrant_is_parent(q1,q2)!=0
end 

function pXest_quadrant_is_parent(::P8estType, q1, q2)
  p8est_quadrant_is_parent(q1,q2)!=0
end 

function pXest_quadrant_is_equal(::P4estType,  q1, q2)
  p4est_quadrant_is_equal(q1, q2)!=0
end 

function pXest_quadrant_is_equal(::P6estType,  q1, q2)
  Gridap.Helpers.@notimplemented
end 

function pXest_quadrant_is_equal(::P8estType,  q1, q2)
  p8est_quadrant_is_equal(q1, q2)!=0
end

function pXest_cell_coords(::P4estType, q, l)
  (q.x,q.y)
end 

function pXest_cell_coords(::P8estType, q, l)
  (q.x,q.y,q.z)
end 

function pXest_cell_coords(::P6estType, q, l)
  (q.x,q.y,l.z)
end

function pXest_num_quadrant_layers(::P4P8estType,q)
  1
end 

function P6EST_COLUMN_GET_RANGE(q)
  f=q.p.piggy3.local_num 
  l=f+q.p.piggy3.which_tree
  (f,l)
end 

function pXest_num_quadrant_layers(::P6estType,q)
  f,l=P6EST_COLUMN_GET_RANGE(q)
  l-f
end 

function pXest_get_layer(::P4P8estType,q,pXest,i)
  q
end 

function p2est_quadrant_array_index(sc_array_object::sc_array_t, it)
  @assert sc_array_object.elem_size == sizeof(p2est_quadrant_t)
  @assert it in 0:sc_array_object.elem_count
  return Ptr{p2est_quadrant_t}(sc_array_object.array + sc_array_object.elem_size*it)
end

function pXest_get_layer(::P6estType,q,pXest,i)
  f,_=P6EST_COLUMN_GET_RANGE(q)
  p2est_quadrant_array_index(pXest.layers[], f+i)[] 
end

function pXest_get_quadrant_and_layer_levels(::P4P8estType,q,l)
  (q.level,)
end 

function pXest_get_quadrant_and_layer_levels(::P6estType,q,l)
  (q.level,l.level)
end 

function pXest_get_quadrant_vertex_coordinates(::P4estType,
                                               connectivity::Ptr{p4est_connectivity_t},
                                               treeid::p4est_topidx_t,
                                               coords,
                                               levels,
                                               corner::Cint,
                                               pvxy::Ptr{Cdouble})
  x,y=coords 
  level,=levels
  p4est_get_quadrant_vertex_coordinates(connectivity,
                                        treeid,
                                        x,
                                        y,
                                        level,
                                        corner,
                                        pvxy)
end

function pXest_get_quadrant_vertex_coordinates(::P6estType,
                                               connectivity::Ptr{p6est_connectivity_t},
                                               treeid::p4est_topidx_t,
                                               coords,
                                               levels,
                                               corner::Cint,
                                               pvxy::Ptr{Cdouble})
  x,y,z=coords 
  qlevel,zlevel=levels
  p6est_get_quadrant_vertex_coordinates(connectivity,
                                        treeid,
                                        x,
                                        y,
                                        z,
                                        qlevel,
                                        zlevel,
                                        corner,
                                        pvxy)
end

function pXest_get_quadrant_vertex_coordinates(::P8estType,
                                               connectivity::Ptr{p8est_connectivity_t},
                                               treeid::p4est_topidx_t,
                                               coords,
                                               levels,
                                               corner::Cint,
                                               pvxy::Ptr{Cdouble})
  x,y,z=coords
  level,=levels 
  p8est_get_quadrant_vertex_coordinates(connectivity,
                                        treeid,
                                        x,
                                        y,
                                        z,
                                        level,
                                        corner,
                                        pvxy)
end


 function _fill_ghost_cells_node_coordinates!(pXest_type::P6estType,
                                             PXEST_CORNERS,
                                             vxy,
                                             pvxy,
                                             node_coordinates,
                                             current,
                                             cell_lids,
                                             ptr_pXest_connectivity,
                                             pXest_ghost)


    function sc_array_p4est_locidx_t_index(sc_array_object::sc_array_t, it)
      @assert sc_array_object.elem_size == sizeof(p4est_locidx_t)
      @assert it in 0:sc_array_object.elem_count
      ptr=Ptr{p4est_locidx_t}(sc_array_object.array + sc_array_object.elem_size*it)
      return unsafe_wrap(Array, ptr, 1)[]
    end

     Dc = num_cell_dims(pXest_type)

     column_ghost = pXest_ghost.column_ghost[]
     ptr_p2est_ghost_quadrants = _unwrap_ghost_quadrants(pXest_type, pXest_ghost)
     ptr_p4est_ghost_quadrants = _unwrap_ghost_quadrants(P4estType(), column_ghost)

     tree_offsets = unsafe_wrap(Array, column_ghost.tree_offsets, pXest_ghost.num_trees+1)
     
     current_ghost_column=0

     # Go over ghost cells
     for i=1:pXest_ghost.num_trees
       for j=tree_offsets[i]:tree_offsets[i+1]-1
          p4est_quadrant = ptr_p4est_ghost_quadrants[j+1]
          k = sc_array_p4est_locidx_t_index(pXest_ghost.column_layer_offsets[],current_ghost_column)
          l = sc_array_p4est_locidx_t_index(pXest_ghost.column_layer_offsets[],current_ghost_column+1)
          for m=k:l-1
            p2est_quadrant = ptr_p2est_ghost_quadrants[m+1]
            coords=pXest_cell_coords(pXest_type,p4est_quadrant,p2est_quadrant)
            levels=pXest_get_quadrant_and_layer_levels(pXest_type,p4est_quadrant,p2est_quadrant)
            for vertex=1:PXEST_CORNERS
                pXest_get_quadrant_vertex_coordinates(pXest_type,
                                                      ptr_pXest_connectivity,
                                                      p4est_topidx_t(i-1),
                                                      coords,
                                                      levels,
                                                      Cint(vertex-1),
                                                      pvxy)
                node_coordinates[cell_lids[current]]=Point{Dc,Float64}(vxy...)
                current=current+1
            end
          end
          current_ghost_column=current_ghost_column+1
       end
     end
 end 

function _fill_ghost_cells_node_coordinates!(pXest_type::P4P8estType,
                                             PXEST_CORNERS,
                                             vxy,
                                             pvxy,
                                             node_coordinates,
                                             current,
                                             cell_lids,
                                             ptr_pXest_connectivity,
                                             pXest_ghost)

     Dc = num_cell_dims(pXest_type)

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
               node_coordinates[cell_lids[current]]=Point{Dc,Float64}(vxy...)
               current=current+1
          end
      end
     end
end 


function generate_node_coordinates(pXest_type::PXestType,
                                   cell_vertex_lids,
                                   nlvertices,
                                   ptr_pXest_connectivity,
                                   ptr_pXest,
                                   ptr_pXest_ghost)

  Dc = num_cell_dims(pXest_type)
  
  PXEST_CORNERS=2^Dc
  pXest_ghost = ptr_pXest_ghost[]
  pXest       = ptr_pXest[]

  dnode_coordinates=map(cell_vertex_lids,nlvertices) do cell_vertex_lids, nl
     node_coordinates=Vector{Point{Dc,Float64}}(undef,nl)
     current=1
     vxy=Vector{Cdouble}(undef,Dc)
     pvxy=pointer(vxy,1)
     cell_lids=cell_vertex_lids.data
     for itree=1:pXest_ghost.num_trees
       tree = pXest_tree_array_index(pXest_type, pXest, itree-1)[]
       # Loop over quadrants/columns in the current tree
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


              @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))] quadrant=$(cell) layer=$(l) coords=$(coords) levels=$(levels) pvxy=$(unsafe_wrap(Array, pvxy, 3)) cell_lids=$(cell_lids[current])"

              node_coordinates[cell_lids[current]]=Point{Dc,Float64}(vxy...)
              current=current+1
            end
          end 
       end
     end
     _fill_ghost_cells_node_coordinates!(pXest_type,
                                         PXEST_CORNERS,
                                         vxy,
                                         pvxy,
                                         node_coordinates,
                                         current,
                                         cell_lids,
                                         ptr_pXest_connectivity,
                                         pXest_ghost)

     node_coordinates
  end
end

function generate_grid_and_topology(pXest_type::PXestType,
                                    cell_vertex_lids,
                                    nlvertices,
                                    node_coordinates)
  Dc=num_cell_dims(pXest_type)
  grid,topology=
  map(cell_vertex_lids,nlvertices,node_coordinates) do cell_vertex_lids, nl, node_coordinates
    polytope= Dc==2 ? QUAD : HEX
    scalar_reffe=Gridap.ReferenceFEs.ReferenceFE(polytope,Gridap.ReferenceFEs.lagrangian,Float64,1)
    cell_types=collect(Fill(1,length(cell_vertex_lids)))
    cell_reffes=[scalar_reffe]
    cell_vertex_lids_gridap=Gridap.Arrays.Table(cell_vertex_lids.data,cell_vertex_lids.ptrs)
    grid = Gridap.Geometry.UnstructuredGrid(node_coordinates,
                                            cell_vertex_lids_gridap,
                                            cell_reffes,
                                            cell_types,
                                            Gridap.Geometry.NonOriented())

    topology = Gridap.Geometry.UnstructuredGridTopology(node_coordinates,
                                      cell_vertex_lids_gridap,
                                      cell_types,
                                      map(Gridap.ReferenceFEs.get_polytope, cell_reffes),
                                      Gridap.Geometry.NonOriented())
    grid,topology
  end |> tuple_of_arrays
  grid,topology
end

function pXest_comm_find_owner(::P4estType,ptr_pXest,itree,quad,guess)
  return p4est_comm_find_owner(ptr_pXest,itree,quad,guess)
end

function pXest_comm_find_owner(::P8estType,ptr_pXest,itree,quad,guess)
  return p8est_comm_find_owner(ptr_pXest,itree,quad,guess)
end

function pXest_compute_migration_control_data(pXest_type::P4P8estType,ptr_pXest_old,ptr_pXest_new)
  pXest_old   = ptr_pXest_old[]
  pXest_new   = ptr_pXest_new[]
  num_trees   = Cint(pXest_old.connectivity[].num_trees)
  my_rank     = pXest_old.mpirank
  ranks_count = Dict{Int,Int}()
  lst_ranks   = Int[]
  old2new     = Vector{Int}(undef,pXest_old.local_num_quadrants)
  current_old_quad_index = 1

  for itree = 0:num_trees-1
    tree = pXest_tree_array_index(pXest_type,pXest_old,itree)[]
    num_quads = Cint(tree.quadrants.elem_count)

    for iquad = 0:num_quads-1
      q = pXest_quadrant_array_index(pXest_type, tree, iquad)
      new_rank = pXest_comm_find_owner(pXest_type,ptr_pXest_new,itree,q,0)
      if (new_rank != my_rank)
        if (!(new_rank+1 in keys(ranks_count)))
          push!(lst_ranks,new_rank+1)
          ranks_count[new_rank+1] = 0
        end
        ranks_count[new_rank+1] += 1
        old2new[current_old_quad_index] = 0
      else
        current_new_quad_index = 1
        new_tree = pXest_tree_array_index(pXest_type,pXest_new,pXest_new.first_local_tree)[]
        for t = pXest_new.first_local_tree:pXest_new.last_local_tree
          new_tree = pXest_tree_array_index(pXest_type,pXest_new,t)[]
          if t == itree
            break
          end
          current_new_quad_index += Cint(new_tree.quadrants.elem_count)
        end
        found = false
        num_quads_new = Cint(new_tree.quadrants.elem_count)
        for iquad_new = 0:num_quads_new-1
          q_new = pXest_quadrant_array_index(pXest_type, new_tree, iquad_new)
          found = pXest_quadrant_is_equal(pXest_type,q,q_new)
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

function pXest_compute_migration_control_data(pXest_type::P6estType,ptr_pXest_old,ptr_pXest_new)
  
  pXest_old = ptr_pXest_old[]
  pXest_new = ptr_pXest_new[]
  
  lst_ranks, columns_snd_lids, columns_old2new=
       pXest_compute_migration_control_data(P4estType(),
                                            pXest_old.columns,
                                            pXest_new.columns)

  ptrs=Vector{Int32}(undef,length(columns_snd_lids.ptrs))
  ptrs.=0
  
  lids=Int32[]
  old2new=Int32[]

  current_col  = 1
  current_old_cell = 1
  current_new_cell = 1
  current_rank = 1 

  num_trees = Cint(pXest_old.columns[].connectivity[].num_trees)
  @assert num_trees == Cint(pXest_new.columns[].connectivity[].num_trees)

  # Go over trees 
  for itree = 0:num_trees-1
    tree = pXest_tree_array_index(pXest_type,pXest_old,itree)[]
    num_quads = Cint(tree.quadrants.elem_count)
    
    # Go over columns of current tree 
    for iquad=0:num_quads-1
      q = pXest_quadrant_array_index(pXest_type,tree,iquad)
      f,l=P6EST_COLUMN_GET_RANGE(q[])
      num_quads_in_column=l-f
      current_new_col=columns_old2new[current_col]
      if (current_new_col==0)
          if ((current_rank+1)!=length(columns_snd_lids.ptrs))
            if (current_col==columns_snd_lids.data[columns_snd_lids.ptrs[current_rank+1]])
              current_rank=current_rank+1
            end 
          end
          ptrs[current_rank+1]+=num_quads_in_column
          for i=0:num_quads_in_column-1
            push!(lids,current_old_cell+i)
            push!(old2new,0)
          end
      else
        # Count how many quads are in previous columns! 
        current_new_cell=1
        col=1
        found=false
        # Go over trees 
        for jtree = 0:num_trees-1
           tree_new = pXest_tree_array_index(pXest_type,pXest_new,jtree)[]
           num_quads_new = Cint(tree_new.quadrants.elem_count)
           # Go over columns of current tree 
           for jquad=0:num_quads_new-1
              if (col==current_new_col)
                found=true 
                break
              end 
              q = pXest_quadrant_array_index(pXest_type,tree_new,jquad)
              f,l=P6EST_COLUMN_GET_RANGE(q[])
              col+=1
              current_new_cell+=l-f
           end
           if (found)
              break
           end
        end
        @assert found
        for i=1:num_quads_in_column
           push!(old2new,current_new_cell)
           current_new_cell+=1
        end
      end
      current_old_cell+=num_quads_in_column
      current_col+=1
    end
  end 
  Gridap.Arrays.length_to_ptrs!(ptrs)
  lst_ranks,PartitionedArrays.JaggedArray(lids,ptrs),old2new
end

function pXest_deflate_quadrants(::P4estType,ptr_pXest,data)
  P4est_wrapper.p4est_deflate_quadrants(ptr_pXest,data)
end

function pXest_deflate_quadrants(::P8estType,ptr_pXest,data)
  P4est_wrapper.p8est_deflate_quadrants(ptr_pXest,data)
end

function pXest_comm_count_pertree(::P4estType,ptr_pXest,pertree)
  p4est_comm_count_pertree(ptr_pXest,pertree)
end

function pXest_comm_count_pertree(::P8estType,ptr_pXest,pertree)
  p8est_comm_count_pertree(ptr_pXest,pertree)
end

function pXest_inflate(::P4estType,
                       comm,
                       ptr_pXest_conn,
                       global_first_quadrant, 
                       pertree,
                       quadrants, 
                       data, 
                       user_pointer)
    P4est_wrapper.p4est_inflate(comm,
                                ptr_pXest_conn,
                                global_first_quadrant, 
                                pertree,
                                quadrants, 
                                data, 
                                user_pointer)
end

function pXest_inflate(::P8estType,
                       comm,
                       ptr_pXest_conn,
                       global_first_quadrant, 
                       pertree,
                       quadrants, 
                       data, 
                       user_pointer)
    P4est_wrapper.p8est_inflate(comm,
                                ptr_pXest_conn,
                                global_first_quadrant, 
                                pertree,
                                quadrants, 
                                data, 
                                user_pointer)
end

function pXest_stride_among_children(::P4P8estType,
                                     ::PXestUniformRefinementRuleType,  
                                     ptr_pXest)
  return 1
end 

function pXest_stride_among_children(::P6estType,
                                     ::PXestVerticalRefinementRuleType,  
                                     ptr_pXest)
  return 1
end

function pXest_stride_among_children(pXest_type::P6estType,
                                     ::PXestHorizontalRefinementRuleType,  
                                     ptr_pXest)
 # Here we are assuming:
 # (1) Each processor has at least one column.
 # (2) The number of layers in each column is the same within and accross processors.
 num_trees = ptr_pXest[].columns[].connectivity[].num_trees
 for itree = 0:num_trees-1
   tree = pXest_tree_array_index(pXest_type, ptr_pXest[], itree)[]
   if tree.quadrants.elem_count>0
     q = pXest_quadrant_array_index(pXest_type, tree, 0)
     f,l=P6EST_COLUMN_GET_RANGE(q[])
     return l-f
   end
 end
 return 0
end

function pXest_uniformly_refine!(::P4estType, ptr_pXest)
  # Refine callbacks
  function refine_fn_2d(::Ptr{p4est_t},which_tree::p4est_topidx_t,quadrant::Ptr{p4est_quadrant_t})
    return Cint(1)
  end
  refine_fn_c = @cfunction($refine_fn_2d,Cint,(Ptr{p4est_t}, p4est_topidx_t, Ptr{p4est_quadrant_t}))
  p4est_refine(ptr_pXest, Cint(0), refine_fn_c, C_NULL)
end

function pXest_uniformly_refine!(::P8estType, ptr_pXest)
  # Refine callbacks
  function refine_fn_3d(::Ptr{p8est_t},which_tree::p4est_topidx_t,quadrant::Ptr{p8est_quadrant_t})
    return Cint(1)
  end
  # C-callable refine callback 3D
  refine_fn_c = @cfunction($refine_fn_3d,Cint,(Ptr{p8est_t}, p4est_topidx_t, Ptr{p8est_quadrant_t}))
  p8est_refine(ptr_pXest, Cint(0), refine_fn_c, C_NULL)
end

function pXest_uniformly_coarsen!(::P4estType, ptr_pXest)
  function coarsen_fn_2d(::Ptr{p4est_t},::p4est_topidx_t,::Ptr{Ptr{p4est_quadrant_t}})
    return Cint(1)
  end
  coarsen_fn_c=@cfunction($coarsen_fn_2d,Cint,(Ptr{p4est_t}, p4est_topidx_t, Ptr{Ptr{p4est_quadrant_t}}))
  p4est_coarsen(ptr_pXest, Cint(0), coarsen_fn_c, C_NULL)
end

function pXest_uniformly_coarsen!(::P8estType, ptr_pXest)
  function coarsen_fn_3d(::Ptr{p8est_t},::p4est_topidx_t,::Ptr{Ptr{p8est_quadrant_t}})
    return Cint(1)
  end
  coarsen_fn_c=@cfunction($coarsen_fn_3d,Cint,(Ptr{p8est_t}, p4est_topidx_t, Ptr{Ptr{p8est_quadrant_t}}))
  p8est_coarsen(ptr_pXest, Cint(0), coarsen_fn_c, C_NULL)
end

function pXest_partition_given!(::P4estType, ptr_pXest, new_num_cells_per_part)
  p4est_partition_given(ptr_pXest, new_num_cells_per_part)
end

function pXest_partition_given!(::P8estType, ptr_pXest, new_num_cells_per_part)
  p8est_partition_given(ptr_pXest, new_num_cells_per_part)
end
