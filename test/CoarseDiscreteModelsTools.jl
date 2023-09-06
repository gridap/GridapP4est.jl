  function setup_model(::Type{Val{3}}, perm)
  #               5 +--------+ 7 
  #                /        /|
  #               /        / |
  #            6 +--------+  |
  #              |        |  |
  #              |  1     |  + 3 
  #              |        | /
  #              |        |/
  #            2 +--------+ 4

  #     6  +--------+ 8 
  #       /        /|
  #      /        / |
  #  11 +--------+  |
  #     |        |  |
  #     |  2     |  + 4
  #     |        | /
  #     |        |/
  #   9 +--------+ 10
    ptr  = [ 1, 9, 17 ]
    if (perm==1)
      data = [ 1,2,3,4,5,6,7,8, 2,9,4,10,6,11,8,12 ]
    elseif (perm==2)
      data = [ 1,2,3,4,5,6,7,8, 10,12,4,8,9,11,2,6 ]
    elseif (perm==3)
      data = [ 1,2,3,4,5,6,7,8, 12,11,8,6,10,9,4,2 ]
    elseif (perm==4) 
      data = [ 1,2,3,4,5,6,7,8, 11,9,6,2,12,10,8,4 ]
    end  
    cell_vertex_lids = Gridap.Arrays.Table(data,ptr)
    node_coordinates = Vector{Point{3,Float64}}(undef,12)
    node_coordinates[1]=Point{3,Float64}(0.0,0.0,0.0)
    node_coordinates[2]=Point{3,Float64}(1.0,0.0,0.0)
    node_coordinates[3]=Point{3,Float64}(0.0,1.0,0.0)
    node_coordinates[4]=Point{3,Float64}(1.0,1.0,0.0)
    node_coordinates[5]=Point{3,Float64}(0.0,0.0,1.0)
    node_coordinates[6]=Point{3,Float64}(1.0,0.0,1.0)
    node_coordinates[7]=Point{3,Float64}(0.0,1.0,1.0)
    node_coordinates[8]=Point{3,Float64}(1.0,1.0,1.0)
    node_coordinates[9]=Point{3,Float64}(2.0,0.0,0.0)
    node_coordinates[10]=Point{3,Float64}(2.0,1.0,0.0)
    node_coordinates[11]=Point{3,Float64}(2.0,0.0,1.0)
    node_coordinates[12]=Point{3,Float64}(2.0,1.0,1.0)

    polytope=HEX
    scalar_reffe=Gridap.ReferenceFEs.ReferenceFE(polytope,Gridap.ReferenceFEs.lagrangian,Float64,1)
    cell_types=collect(Fill(1,length(cell_vertex_lids)))
    cell_reffes=[scalar_reffe]
    grid = Gridap.Geometry.UnstructuredGrid(node_coordinates,
                                            cell_vertex_lids,
                                            cell_reffes,
                                            cell_types,
                                            Gridap.Geometry.NonOriented())
    m=Gridap.Geometry.UnstructuredDiscreteModel(grid)
    labels = get_face_labeling(m)
    labels.d_to_dface_to_entity[1].=2
    if (perm==1 || perm==2)
      labels.d_to_dface_to_entity[2].=2
      labels.d_to_dface_to_entity[3].=[2,2,2,2,2,1,2,2,2,2,2]
    elseif (perm==3 || perm==4)
      labels.d_to_dface_to_entity[2].=2
      labels.d_to_dface_to_entity[3].=[2,2,2,2,2,1,2,2,2,2,2]
    end
    labels.d_to_dface_to_entity[4].=1  
    add_tag!(labels,"boundary",[2])
    add_tag!(labels,"interior",[1])
    m
  end

  function setup_model(::Type{Val{2}}, perm)
    @assert perm ∈ (1,2,3,4)
    #
    #  3-------4-------6
    #  |       |       |
    #  |       |       |
    #  |       |       |
    #  1-------2-------5
    #
        ptr  = [ 1, 5, 9 ]
        if (perm==1)
          data = [ 1,2,3,4, 2,5,4,6 ]
        elseif (perm==2)
          data = [ 1,2,3,4, 6,4,5,2 ]
        elseif (perm==3)
          data = [ 4,3,2,1, 2,5,4,6 ]
        elseif (perm==4) 
          data = [ 4,3,2,1, 6,4,5,2 ]
        end  
        cell_vertex_lids = Gridap.Arrays.Table(data,ptr)
        node_coordinates = Vector{Point{2,Float64}}(undef,6)
        node_coordinates[1]=Point{2,Float64}(0.0,0.0)
        node_coordinates[2]=Point{2,Float64}(1.0,0.0)
        node_coordinates[3]=Point{2,Float64}(0.0,1.0)
        node_coordinates[4]=Point{2,Float64}(1.0,1.0)
        node_coordinates[5]=Point{2,Float64}(2.0,0.0)
        node_coordinates[6]=Point{2,Float64}(2.0,1.0)

        polytope=QUAD
        scalar_reffe=Gridap.ReferenceFEs.ReferenceFE(polytope,Gridap.ReferenceFEs.lagrangian,Float64,1)
        cell_types=collect(Fill(1,length(cell_vertex_lids)))
        cell_reffes=[scalar_reffe]
        grid = Gridap.Geometry.UnstructuredGrid(node_coordinates,
                                                cell_vertex_lids,
                                                cell_reffes,
                                                cell_types,
                                                Gridap.Geometry.NonOriented())
        m=Gridap.Geometry.UnstructuredDiscreteModel(grid)
        labels = get_face_labeling(m)
        labels.d_to_dface_to_entity[1].=2
        if (perm==1 || perm==2)
          labels.d_to_dface_to_entity[2].=[2,2,2,1,2,2,2]
        elseif (perm==3 || perm==4)
          labels.d_to_dface_to_entity[2].=[2,2,1,2,2,2,2] 
        end
        labels.d_to_dface_to_entity[3].=1
        add_tag!(labels,"boundary",[2])
        add_tag!(labels,"interior",[1])
        m
  end