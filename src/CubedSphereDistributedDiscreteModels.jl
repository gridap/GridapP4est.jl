using Gridap
using GridapDistributed
using MPI
using p4est_wrapper
using FillArrays

const ref_panel_coordinates = [ -1.0, -1.0, 1.0, -1.0, -1.0, 1.0, 1.0, 1.0 ]

function set_panel_vertices_αβ_coordinates!( pconn :: Ptr{p4est_wrapper.p4est_connectivity_t},
                                             coarse_discrete_model :: DiscreteModel{2,2},
                                             panel )
  @assert num_cells(coarse_discrete_model) == 6
  @assert panel ≤ num_cells(coarse_discrete_model)
  @assert panel ≥ 1
  trian=Triangulation(coarse_discrete_model)
  cell_vertices=Gridap.Geometry.get_cell_node_ids(trian)
  #println(cell_vertices)
  cell_vertices_panel=cell_vertices[panel]
  conn=pconn[]
  vertices=unsafe_wrap(Array,
                       conn.vertices,
                       length(Gridap.Geometry.get_node_coordinates(coarse_discrete_model))*3)
  for (l,g) in enumerate(cell_vertices_panel)
     # println("XXX $(l) $(g) $(lref)")
     vertices[(g-1)*3+1]=ref_panel_coordinates[(l-1)*2+1]
     vertices[(g-1)*3+2]=ref_panel_coordinates[(l-1)*2+2]
  end
end

function setup_cubed_sphere_coarse_discrete_model()
    # 6 panels (cells), 4 corners (vertices) each panel
    ptr  = [ 1, 5, 9, 13, 17, 21, 25 ]
    data = [ 1,2,3,4, 2,5,4,6, 7,8,5,6, 1,3,7,8, 3,4,8,6, 1,7,2,5  ]
    cell_vertex_lids = Gridap.Arrays.Table(data,ptr)
    node_coordinates = Vector{Point{2,Float64}}(undef,8)
    for i in 1:length(node_coordinates)
      node_coordinates[i]=Point{2,Float64}(0.0,0.0)
    end

    polytope=QUAD
    scalar_reffe=Gridap.ReferenceFEs.ReferenceFE(polytope,Gridap.ReferenceFEs.lagrangian,Float64,1)
    cell_types=collect(Fill(1,length(cell_vertex_lids)))
    cell_reffes=[scalar_reffe]
    grid = Gridap.Geometry.UnstructuredGrid(node_coordinates,
                                            cell_vertex_lids,
                                            cell_reffes,
                                            cell_types,
                                            Gridap.Geometry.NonOriented())
    Gridap.Geometry.UnstructuredDiscreteModel(grid)
end

function generate_ptr(ncells)
  nvertices=4
  ptr  = Vector{Int}(undef,ncells+1)
  ptr[1]=1
  for i=1:ncells
    ptr[i+1]=ptr[i]+nvertices
  end
  ptr
end

function generate_cell_coordinates_and_panels(comm,
                                   coarse_discrete_model,
                                   ptr_pXest_connectivity,
                                   ptr_pXest,
                                   ptr_pXest_ghost)

  Dc=2
  PXEST_CORNERS=4
  pXest_ghost = ptr_pXest_ghost[]
  pXest = ptr_pXest[]

  # Obtain ghost quadrants
  ptr_ghost_quadrants = Ptr{p4est_wrapper.p4est_quadrant_t}(pXest_ghost.ghosts.array)

  tree_offsets = unsafe_wrap(Array, pXest_ghost.tree_offsets, pXest_ghost.num_trees+1)
  dcell_coordinates_and_panels=DistributedData(comm) do part
     ncells=pXest.local_num_quadrants+pXest_ghost.ghosts.elem_count
     panels = Vector{Int}(undef,ncells)
     data = Vector{Point{Dc,Float64}}(undef,ncells*PXEST_CORNERS)
     ptr  = generate_ptr(ncells)
     current=1
     current_cell=1
     vxy=Vector{Cdouble}(undef,Dc)
     pvxy=pointer(vxy,1)
     for itree=1:pXest_ghost.num_trees
       tree = p4est_tree_array_index(pXest.trees, itree-1)[]
       if tree.quadrants.elem_count > 0
          set_panel_vertices_αβ_coordinates!( ptr_pXest_connectivity, coarse_discrete_model, itree)
       end
       for cell=1:tree.quadrants.elem_count
          panels[current_cell]=itree
          quadrant=p4est_quadrant_array_index(tree.quadrants, cell-1)[]
          for vertex=1:PXEST_CORNERS
            GridapDistributed.p4est_get_quadrant_vertex_coordinates(ptr_pXest_connectivity,
                                                  p4est_topidx_t(itree-1),
                                                  quadrant.x,
                                                  quadrant.y,
                                                  quadrant.level,
                                                  Cint(vertex-1),
                                                  pvxy)
            data[current]=Point{Dc,Float64}(vxy...)
            current=current+1
          end
          current_cell=current_cell+1
       end
     end

     # Go over ghost cells
     for i=1:pXest_ghost.num_trees
      if tree_offsets[i+1]-tree_offsets[i] > 0
        set_panel_vertices_αβ_coordinates!( ptr_pXest_connectivity, coarse_discrete_model, i)
      end
      for j=tree_offsets[i]:tree_offsets[i+1]-1
          panels[current_cell]=i
          quadrant = ptr_ghost_quadrants[j+1]
          for vertex=1:PXEST_CORNERS
            GridapDistributed.p4est_get_quadrant_vertex_coordinates(ptr_pXest_connectivity,
                                                     p4est_topidx_t(i-1),
                                                     quadrant.x,
                                                     quadrant.y,
                                                     quadrant.level,
                                                     Cint(vertex-1),
                                                     pvxy)

          #  if (MPI.Comm_rank(comm.comm)==0)
          #     println(vxy)
          #  end
          data[current]=Point{Dc,Float64}(vxy...)
          current=current+1
         end
         current_cell=current_cell+1
       end
     end
     Gridap.Arrays.Table(data,ptr), panels
  end
end

function generate_αβgrid_geo(cell_coordinates_and_panels)
  function map_panel_αβ_2_xyz(αβ,panel)
    α,β=αβ
    if panel==1
      x=Point(1.0,α,β)
    elseif panel==2
      x=Point(-α,1.0,β)
    elseif panel==3
      x=Point(-1.0,β,α)
    elseif panel==4
      x=Point(-β,-1.0,α)
    elseif panel==5
      x=Point(-β,α,1.0)
    elseif panel==6
      x=Point(-α,β,-1.0)
    end
    x
  end

  DistributedData(cell_coordinates_and_panels) do part, (cell_coordinates,panels)
     ptr  = generate_ptr(length(cell_coordinates))
     data = collect(1:length(cell_coordinates)*4)
     cell_vertex_lids=Gridap.Arrays.Table(data,ptr)
     node_coordinates=Vector{Point{3,Float64}}(undef,length(cell_coordinates)*4)

     cache=array_cache(cell_coordinates)
     current=1
     for cell=1:length(cell_coordinates)
        current_cell_coordinates=getindex!(cache,cell_coordinates,cell)
        for i=1:length(current_cell_coordinates)
          node_coordinates[current]=map_panel_αβ_2_xyz(current_cell_coordinates[i],panels[cell])
          current=current+1
        end
     end

     polytope=QUAD
     scalar_reffe=Gridap.ReferenceFEs.ReferenceFE(polytope,Gridap.ReferenceFEs.lagrangian,Float64,1)
     cell_types=collect(Fill(1,length(cell_vertex_lids)))
     cell_reffes=[scalar_reffe]
     grid=Gridap.Geometry.UnstructuredGrid(node_coordinates,
                                      cell_vertex_lids,
                                      cell_reffes,
                                      cell_types,
                                      Gridap.Geometry.NonOriented())
     grid
  end
end

function generate_αβgrid_top(cell_vertex_lids_nlvertices)
  DistributedData(cell_vertex_lids_nlvertices) do part, (cell_vertex_lids,nlvector)
     node_coordinates=Vector{Point{2,Float64}}(undef,nlvector)
     polytope=QUAD
     scalar_reffe=Gridap.ReferenceFEs.ReferenceFE(polytope,Gridap.ReferenceFEs.lagrangian,Float64,1)
     cell_types=collect(Fill(1,length(cell_vertex_lids)))
     cell_reffes=[scalar_reffe]
    #  if (part==2)
    #   println(cell_vertex_lids)
    #  end
     grid=Gridap.Geometry.UnstructuredGrid(node_coordinates,
                                      cell_vertex_lids,
                                      cell_reffes,
                                      cell_types,
                                      Gridap.Geometry.NonOriented())
     grid
  end
end

function CubedSphereDiscreteModel(comm::Communicator,num_uniform_refinements::Int)

  Dc=2

  coarse_discrete_model=setup_cubed_sphere_coarse_discrete_model()

  ptr_pXest_connectivity=GridapDistributed.setup_pXest_connectivity(coarse_discrete_model)

  # Create a new forest
  ptr_pXest = GridapDistributed.setup_pXest(Val{Dc},comm,ptr_pXest_connectivity,num_uniform_refinements)

  # Build the ghost layer
  ptr_pXest_ghost=GridapDistributed.setup_pXest_ghost(Val{Dc},ptr_pXest)

  cellindices = GridapDistributed.setup_cell_indexset(Val{Dc},comm,ptr_pXest,ptr_pXest_ghost)

  ptr_pXest_lnodes=GridapDistributed.setup_pXest_lnodes(Val{Dc}, ptr_pXest, ptr_pXest_ghost)

  cell_vertex_gids=GridapDistributed.generate_cell_vertex_gids(ptr_pXest_lnodes,cellindices)

  cell_vertex_lids_nlvertices=GridapDistributed.generate_cell_vertex_lids_nlvertices(cell_vertex_gids)

  cell_coordinates_and_panels=generate_cell_coordinates_and_panels(comm,
                                             coarse_discrete_model,
                                             ptr_pXest_connectivity,
                                             ptr_pXest,
                                             ptr_pXest_ghost)

  αβgrid_geo=generate_αβgrid_geo(cell_coordinates_and_panels)
  αβgrid_top=generate_αβgrid_top(cell_vertex_lids_nlvertices)

  # do_on_parts(comm, cell_coordinates) do part, cell_coords
  #   if part==1
  #     println(cell_coords[9:12])
  #   end
  # end

  ddiscretemodel=
  DistributedData(comm,αβgrid_geo,αβgrid_top) do part, αβgrid_geo, αβgrid_top
    AnalyticalEquiAngularMapCubedSphereDiscreteModel(αβgrid_geo, αβgrid_top)
  end

  # Write forest to VTK file
  #p4est_vtk_write_file(unitsquare_forest, C_NULL, "my_step")

  # Destroy lnodes
  p4est_ghost_destroy(ptr_pXest_ghost)
  # Destroy the forest
  p4est_destroy(ptr_pXest)
  # Destroy the connectivity
  p4est_connectivity_destroy(ptr_pXest_connectivity)

  GridapDistributed.DistributedDiscreteModel(ddiscretemodel,cellindices)
end


struct MapCubeToSphere{T} <: Function
  radius::T
end

function (map::MapCubeToSphere{T})(xyz) where T
  x,y,z = xyz
  xₛ = x*sqrt(1.0-y^2/2-z^2/2+y^2*z^2/(3.0))
  yₛ = y*sqrt(1.0-z^2/2-x^2/2+x^2*z^2/(3.0))
  zₛ = z*sqrt(1.0-x^2/2-y^2/2+x^2*y^2/(3.0))
  map.radius*Point(xₛ,yₛ,zₛ)
end



struct AnalyticalEquiAngularMapCubedSphereTriangulation{T} <: Triangulation{2,3}
  cell_map::T
  αβgrid::Gridap.Geometry.UnstructuredGrid{2,3}
end

# Triangulation API

# Delegating to the underlying face Triangulation

Gridap.Geometry.get_cell_coordinates(trian::AnalyticalEquiAngularMapCubedSphereTriangulation) = Gridap.Geometry.get_cell_coordinates(trian.αβgrid)

Gridap.Geometry.get_reffes(trian::AnalyticalEquiAngularMapCubedSphereTriangulation) = Gridap.Geometry.get_reffes(trian.αβgrid)

Gridap.Geometry.get_cell_type(trian::AnalyticalEquiAngularMapCubedSphereTriangulation) = Gridap.Geometry.get_cell_type(trian.αβgrid)

Gridap.Geometry.get_node_coordinates(trian::AnalyticalEquiAngularMapCubedSphereTriangulation) = Gridap.Helpers.get_node_coordinates(trian.αβgrid)

Gridap.Geometry.get_cell_node_ids(trian::AnalyticalEquiAngularMapCubedSphereTriangulation) = Gridap.Helpers.get_cell_node_ids(trian.αβgrid)

Gridap.Geometry.get_cell_map(trian::AnalyticalEquiAngularMapCubedSphereTriangulation) = trian.cell_map

# Genuine methods

Gridap.Geometry.TriangulationStyle(::Type{<:AnalyticalEquiAngularMapCubedSphereTriangulation}) = SubTriangulation()

Gridap.Geometry.get_background_triangulation(trian::AnalyticalEquiAngularMapCubedSphereTriangulation) =
      Gridap.Geometry.get_background_triangulation(trian.αβgrid)

Gridap.Geometry.get_cell_to_bgcell(trian::AnalyticalEquiAngularMapCubedSphereTriangulation) = get_cell_to_bgcell(trian.αβgrid)

function Gridap.Geometry.get_cell_to_bgcell(
  trian_in::AnalyticalEquiAngularMapCubedSphereTriangulation,
  trian_out::AnalyticalEquiAngularMapCubedSphereTriangulation)
  Gridap.Helpers.@notimplemented
end

function Gridap.Geometry.is_included(
  trian_in::AnalyticalEquiAngularMapCubedSphereTriangulation,
  trian_out::AnalyticalEquiAngularMapCubedSphereTriangulation)
  Gridap.Helpers.@notimplemented
end

function Gridap.Geometry.get_cell_ref_map(trian::AnalyticalEquiAngularMapCubedSphereTriangulation)
  get_cell_ref_map(trian.btrian)
end

function Gridap.Geometry.get_cell_ref_map(
  trian_in::AnalyticalEquiAngularMapCubedSphereTriangulation,
  trian_out::AnalyticalEquiAngularMapCubedSphereTriangulation)
  Gridap.Helpers.@notimplemented
end


struct AnalyticalEquiAngularMapCubedSphereDiscreteModel{T,B,C} <: Gridap.Geometry.DiscreteModel{2,3}
  cell_map::T
  cubed_sphere_model::B
  trian::C
  function AnalyticalEquiAngularMapCubedSphereDiscreteModel(
      αβgrid_geo::Gridap.Geometry.UnstructuredGrid{2,3},
      αβgrid_top::Gridap.Geometry.UnstructuredGrid{2,2};
      radius=1)
    m1=Fill(Gridap.Fields.GenericField(MapCubeToSphere(radius)),num_cells(αβgrid_geo))
    m2=get_cell_map(αβgrid_geo)
    m=lazy_map(∘,m1,m2)

    cubed_mesh_model=Gridap.Geometry.UnstructuredDiscreteModel(αβgrid_top)

    # Wrap up BoundaryTriangulation
    trian=AnalyticalEquiAngularMapCubedSphereTriangulation(m,αβgrid_geo)

    # Build output object
    T=typeof(m)
    B=typeof(cubed_mesh_model)
    C=typeof(trian)
    GC.gc()
    new{T,B,C}(m,cubed_mesh_model,trian)
  end
end

Gridap.Geometry.get_cell_map(model::AnalyticalEquiAngularMapCubedSphereDiscreteModel) = model.cell_map
Gridap.Geometry.get_grid(model::AnalyticalEquiAngularMapCubedSphereDiscreteModel) = Gridap.Geometry.get_grid(model.cubed_sphere_model)
Gridap.Geometry.get_grid_topology(model::AnalyticalEquiAngularMapCubedSphereDiscreteModel) = Gridap.Geometry.get_grid_topology(model.cubed_sphere_model)
Gridap.Geometry.get_face_labeling(model::AnalyticalEquiAngularMapCubedSphereDiscreteModel) = Gridap.Geometry.get_face_labeling(model.cubed_sphere_model)
Gridap.Geometry.get_triangulation(a::AnalyticalEquiAngularMapCubedSphereDiscreteModel) = a.trian
Gridap.Geometry.Triangulation(a::AnalyticalEquiAngularMapCubedSphereDiscreteModel) = a.trian

function Gridap.CellData.get_normal_vector(model::AnalyticalEquiAngularMapCubedSphereDiscreteModel)
    cell_normal = Gridap.Geometry.get_facet_normal(model)
    Gridap.CellData.GenericCellField(cell_normal,Triangulation(model),ReferenceDomain())
end

function _unit_outward_normal(v::Gridap.Fields.MultiValue{Tuple{2,3}})
  n1 = v[1,2]*v[2,3] - v[1,3]*v[2,2]
  n2 = v[1,3]*v[2,1] - v[1,1]*v[2,3]
  n3 = v[1,1]*v[2,2] - v[1,2]*v[2,1]
  n = VectorValue(n1,n2,n3)
  n/norm(n)
end


function Gridap.Geometry.get_facet_normal(model::AnalyticalEquiAngularMapCubedSphereDiscreteModel)
  # Get the Jacobian of the cubed sphere mesh
  map   = get_cell_map(model)
  Jt    = lazy_map(∇,map)
  p=lazy_map(Operation(_unit_outward_normal),Jt)
  p
end

function run(comm)
  num_uniform_refinements=4
  model=CubedSphereDiscreteModel(comm,num_uniform_refinements)
  do_on_parts(model) do part, (model,gids)
    writevtk(Triangulation(model),
             "part_$(part)",
             cellfields=["n"=>Gridap.CellData.get_normal_vector(model),
                         "lid_to_owner"=>gids.lid_to_owner])
  end
end

MPIPETScCommunicator() do comm
  run(comm)
end
