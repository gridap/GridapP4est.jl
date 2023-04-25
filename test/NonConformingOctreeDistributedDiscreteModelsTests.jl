#module NonConformingOctreeDistributedDiscreteModelsTests
  using P4est_wrapper
  using GridapP4est
  using Gridap
  using PartitionedArrays
  using GridapDistributed
  using MPI
  using Gridap.FESpaces
  using FillArrays
  # Generate a local numbering of vertices that includes hanging vertices 
  # Generate a local numbering of faces out of the one generated by vertices (automatic? to confirm)

  # Establish the correspondence among local numbering of vertices and p4est global numbering 
  # Establish the correspondence among local numbering of faces and p4est global numbering 

  # Generate a global numbering of (regular,hanging) vertices?
  # Generate a global numbering of (regular,hanging) faces?

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
      @assert false
    elseif (perm==3)
      @assert false
    elseif (perm==4) 
      @assert false
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
        @assert false 
    end  
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
        add_tag!(labels,"boundary",[2])
        add_tag!(labels,"interior",[1])
        m
  end


  ## Better to use a C-enum. But I did not use it in order to keep the Julia
  ## version of this C example as simple as possible
  const nothing_flag = Cint(0)
  const refine_flag = Cint(1)

  ## Refine those cells with even identifier    (0,2,4,6,8,...)
  ## Leave untouched cells with odd identifier  (1,3,5,7,9,...)
  function allocate_and_set_refinement_and_coarsening_flags(forest_ptr::Ptr{p4est_t})
    forest = forest_ptr[]
    tree = p4est_tree_array_index(forest.trees, 0)[]
    return [i != 1 ? nothing_flag : refine_flag for i = 1:tree.quadrants.elem_count]
  end

  MPI.Init()
  parts = get_part_ids(MPIBackend(), 1)
  # run(parts, (1, 1))
  # MPI.Finalize()

  function test(TVDc::Type{Val{Dc}}, perm, order) where Dc
    # This is for debuging
    coarse_model = setup_model(TVDc,perm)
    model = OctreeDistributedDiscreteModel(parts, coarse_model, 0)

    ref_coarse_flags=map_parts(parts) do _
      [refine_flag,nothing_flag]
      #allocate_and_set_refinement_and_coarsening_flags(model.ptr_pXest)
    end 
    dmodel,non_conforming_glue=refine(model,ref_coarse_flags)

    println(non_conforming_glue)

    p8est_vtk_write_file(dmodel.ptr_pXest, C_NULL, string("adapted_forest"))

    # FE Spaces
    reffe = ReferenceFE(lagrangian,Float64,order)
    V = TestFESpace(dmodel,reffe,dirichlet_tags="boundary")
    U = TrialFESpace(V)   

    function _h_refined_reffe(reffe::Tuple{<:Lagrangian,Any,Any})
      (reffe[1],(reffe[2][1],2*reffe[2][2]),reffe[3])
    end

    cell_polytope = Dc==2 ? QUAD : HEX


    basis, reffe_args,reffe_kwargs = reffe
    cell_reffe = ReferenceFE(cell_polytope, basis,reffe_args...;reffe_kwargs...)
    reffe_cell = cell_reffe

    h_refined_reffe = _h_refined_reffe(reffe)
    basis, reffe_args,reffe_kwargs = h_refined_reffe
    cell_reffe_h_refined = ReferenceFE(cell_polytope,basis,reffe_args...;reffe_kwargs...)
    reffe_cell_h_refined = cell_reffe_h_refined

    dof_basis_h_refined = Gridap.CellData.get_dof_basis(reffe_cell_h_refined)

    coarse_shape_funs=Gridap.ReferenceFEs.get_shapefuns(reffe_cell)
    ref_constraints=evaluate(dof_basis_h_refined,coarse_shape_funs)


    # To-think: might this info go to the glue? 
    # If it is required in different scenarios, I would say it may make sense
    function _generate_hanging_faces_to_cell_and_lface(num_regular_faces, 
                                                      num_hanging_faces, 
                                                      gridap_cell_faces)
      # Locate for each hanging vertex a cell to which it belongs 
      # and local position within that cell 
      hanging_faces_to_cell = Vector{Int}(undef, num_hanging_faces) 
      hanging_faces_to_lface = Vector{Int}(undef, num_hanging_faces)
      for cell=1:length(gridap_cell_faces)
        s=gridap_cell_faces.ptrs[cell]
        e=gridap_cell_faces.ptrs[cell+1]
        l=e-s
        for j=1:l
          fid=gridap_cell_faces.data[s+j-1]
          if fid>num_regular_faces
            fid_hanging=fid-num_regular_faces
            hanging_faces_to_cell[fid_hanging]=cell
            hanging_faces_to_lface[fid_hanging]=j
          end 
        end 
      end
      hanging_faces_to_cell, hanging_faces_to_lface
    end


    function _generate_hanging_faces_owner_face_dofs(num_hanging_faces, 
                                                    face_dofs,
                                                    hanging_faces_glue,
                                                    cell_dof_ids)
      
      cache = array_cache(cell_dof_ids)
      ptrs=Vector{Int}(undef, num_hanging_faces+1)
      ptrs[1]=1
      for fid_hanging=1:num_hanging_faces
        glue=hanging_faces_glue[fid_hanging]
        ocell_lface=glue[2]
        ptrs[fid_hanging+1] = ptrs[fid_hanging] + length(face_dofs[ocell_lface])
      end 
      data_owner_face_dofs=Vector{Int}(undef, ptrs[num_hanging_faces+1]-1)
      for fid_hanging=1:num_hanging_faces
        glue=hanging_faces_glue[fid_hanging]
        ocell=glue[1]
        ocell_lface=glue[2]
        s=ptrs[fid_hanging]
        e=ptrs[fid_hanging+1]-1
        current_cell_dof_ids = getindex!(cache, cell_dof_ids, ocell)
        for (j,ldof) in enumerate(face_dofs[ocell_lface])
          data_owner_face_dofs[s+j-1]=current_cell_dof_ids[ldof]
        end
      end 
      Gridap.Arrays.Table(data_owner_face_dofs, ptrs)
    end

    function face_dim(::Type{Val{Dc}}, face_lid) where Dc
      num_vertices = GridapP4est.num_cell_vertices(Val{Dc})
      num_edges    = GridapP4est.num_cell_edges(Val{Dc})
      num_faces    = GridapP4est.num_cell_faces(Val{Dc})
      if (face_lid <= num_vertices)
        return 0
      elseif (face_lid <= num_vertices+num_edges)
        return 1
      elseif (face_lid <= num_vertices+num_edges+num_faces)
        return Dc-1
      end  
    end

    function face_lid_within_dim(::Type{Val{Dc}}, face_lid) where Dc
      num_vertices = GridapP4est.num_cell_vertices(Val{Dc})
      num_edges    = GridapP4est.num_cell_edges(Val{Dc})
      num_faces    = GridapP4est.num_cell_faces(Val{Dc})
      if (face_lid<=num_vertices)
        return face_lid
      elseif (face_lid <= num_vertices+num_edges)
        return face_lid-num_vertices
      elseif (face_lid <= num_vertices+num_edges+num_faces)
        return face_lid-num_vertices-num_edges
      end  
    end

    function _generate_constraints!(Df,
                                    Dc,
                                    cell_faces,
                                    num_hanging_faces,
                                    hanging_faces_to_cell,
                                    hanging_faces_to_lface,
                                    hanging_faces_owner_face_dofs,
                                    hanging_faces_glue,
                                    face_subface_ldof_to_cell_ldof,
                                    face_dofs,
                                    face_own_dofs,
                                    subface_own_dofs,
                                    cell_dof_ids,
                                    node_permutations,
                                    owner_faces_pindex,
                                    owner_faces_lids,
                                    ref_constraints,
                                    sDOF_to_dof,
                                    sDOF_to_dofs,
                                    sDOF_to_coeffs)

      @assert Dc==2 || Dc==3
      @assert 0 ≤ Df < Dc

      num_vertices = GridapP4est.num_cell_vertices(Val{Dc})
      num_edges    = GridapP4est.num_cell_edges(Val{Dc})
      num_faces    = GridapP4est.num_cell_faces(Val{Dc})
      
      offset=0
      if (Df ≥ 1)
        offset+=num_vertices
        if (Df == 2)
          offset+=num_edges
        end 
      end 

      cache_dof_ids=array_cache(cell_dof_ids)
      first_cache_cell_faces=array_cache(first(cell_faces))
      cache_cell_faces=Vector{typeof(first_cache_cell_faces)}(undef, length(cell_faces))
      for i=1:length(cell_faces)
        cache_cell_faces[i]=array_cache(cell_faces[i])
      end

      for fid_hanging=1:num_hanging_faces
        cell=hanging_faces_to_cell[fid_hanging]
        current_cell_dof_ids=getindex!(cache_dof_ids,cell_dof_ids,cell)
        lface=hanging_faces_to_lface[fid_hanging]
        ocell,ocell_lface,subface=hanging_faces_glue[fid_hanging]
        ocell_lface_within_dim = face_lid_within_dim(Val{Dc}, ocell_lface)
        oface_dim = face_dim(Val{Dc}, ocell_lface)
        
        if (Df==0) # Am I a vertex?
          hanging_lvertex_within_first_subface = 2^oface_dim
          cur_subface_own_dofs=subface_own_dofs[oface_dim][hanging_lvertex_within_first_subface]
        elseif (Df==Dc-1) # Am I a face?
          @assert oface_dim == Dc-1
          cur_subface_own_dofs=subface_own_dofs[oface_dim][end]
        end

        
        oface=getindex!(cache_cell_faces[oface_dim+1],cell_faces[oface_dim+1],ocell)[ocell_lface_within_dim]
        oface_lid,_=owner_faces_lids[oface_dim][oface]
        pindex=owner_faces_pindex[oface_dim][oface_lid]
        for ((ldof,dof),ldof_subface) in zip(enumerate(face_own_dofs[offset+lface]),cur_subface_own_dofs)
          push!(sDOF_to_dof,current_cell_dof_ids[dof])
          push!(sDOF_to_dofs,hanging_faces_owner_face_dofs[fid_hanging])
          coeffs=Vector{Float64}(undef,length(hanging_faces_owner_face_dofs[fid_hanging]))
          # Go over dofs of ocell_lface
          for (ifdof,icdof) in enumerate(face_dofs[ocell_lface])
            pifdof=node_permutations[oface_dim][pindex][ifdof]
            println("XXXX: $(ifdof) $(pifdof)")
            ldof_coarse=face_dofs[ocell_lface][pifdof]
            coeffs[ifdof]=
            ref_constraints[face_subface_ldof_to_cell_ldof[oface_dim][ocell_lface_within_dim][subface][ldof_subface],ldof_coarse]
          end 
          push!(sDOF_to_coeffs,coeffs) 
        end 
      end
    end


    # count how many different owner faces
    # for each owner face 
    #    track the global IDs of its face vertices from the perspective of the subfaces
    # for each owner face 
    #    compute permutation id
    function _compute_owner_faces_pindex_and_lids(Df,
                                        Dc,
                                        num_hanging_faces, 
                                        hanging_faces_glue,
                                        hanging_faces_to_cell, 
                                        hanging_faces_to_lface,
                                        cell_vertices,
                                        cell_faces,
                                        lface_to_cvertices,
                                        pindex_to_cfvertex_to_fvertex)
      num_owner_faces=0
      owner_faces_lids=Dict{Int,Tuple{Int,Int,Int}}()
      for fid_hanging=1:num_hanging_faces
        ocell,ocell_lface,_=hanging_faces_glue[fid_hanging]
        ocell_dim = face_dim(Val{Dc}, ocell_lface)
        if (ocell_dim==Df)
          ocell_lface_within_dim = face_lid_within_dim(Val{Dc}, ocell_lface)
          owner_face=cell_faces[ocell][ocell_lface_within_dim]
          if !(haskey(owner_faces_lids,owner_face))
            num_owner_faces+=1
            owner_faces_lids[owner_face]=(num_owner_faces,ocell,ocell_lface)
          end
        end 
      end

      num_face_vertices=length(first(lface_to_cvertices))
      owner_face_vertex_ids=Vector{Int}(undef,num_face_vertices*num_owner_faces)

      for fid_hanging=1:num_hanging_faces
        ocell,ocell_lface,subface=hanging_faces_glue[fid_hanging]
        ocell_dim = face_dim(Val{Dc}, ocell_lface)
        if (ocell_dim==Df)
          cell=hanging_faces_to_cell[fid_hanging]   
          lface=hanging_faces_to_lface[fid_hanging]
          cvertex=lface_to_cvertices[lface][subface]
          vertex=cell_vertices[cell][cvertex]
          ocell_lface_within_dim = face_lid_within_dim(Val{Dc}, ocell_lface)
          owner_face=cell_faces[ocell][ocell_lface_within_dim]
          owner_face_lid,_=owner_faces_lids[owner_face]
          owner_face_vertex_ids[(owner_face_lid-1)*num_face_vertices+subface]=vertex
        end    
      end

      owner_faces_pindex=Vector{Int}(undef,num_owner_faces)
      for owner_face in keys(owner_faces_lids)
        (owner_face_lid,ocell,ocell_lface)=owner_faces_lids[owner_face]
        ocell_lface_within_dim = face_lid_within_dim(Val{Dc}, ocell_lface)
        # Compute permutation id by comparing 
        #  1. cell_vertices[ocell][ocell_lface]
        #  2. owner_face_vertex_ids 
        pindexfound = false
        cfvertex_to_cvertex = lface_to_cvertices[ocell_lface_within_dim]
        for (pindex, cfvertex_to_fvertex) in enumerate(pindex_to_cfvertex_to_fvertex)
          found = true
          for (cfvertex,fvertex) in enumerate(cfvertex_to_fvertex)
            vertex1 = owner_face_vertex_ids[(owner_face_lid-1)*num_face_vertices+fvertex]
            cvertex = cfvertex_to_cvertex[cfvertex]
            vertex2 = cell_vertices[ocell][cvertex]
            if vertex1 != vertex2
              found = false
              break
            end
          end
          if found
            owner_faces_pindex[owner_face_lid] = pindex
            pindexfound = true
            break
          end
        end
        @assert pindexfound "Valid pindex not found"
      end
      owner_faces_pindex, owner_faces_lids
    end 

    function generate_constraints(dmodel::GridapDistributed.DistributedDiscreteModel{Dc}, 
                                  V,
                                  reffe,
                                  non_conforming_glue,
                                  ref_constraints, 
                                  face_subface_ldof_to_cell_ldof) where Dc
      num_regular_faces,
      num_hanging_faces,
      gridap_cell_faces,
      hanging_faces_glue = non_conforming_glue

      gridap_cell_faces = map_parts(gridap_cell_faces...) do gridap_cell_faces...
        gridap_cell_faces
      end 
      num_regular_faces = map_parts(num_regular_faces...) do num_regular_faces...
        num_regular_faces
      end
      num_hanging_faces = map_parts(num_hanging_faces...) do num_hanging_faces...
        num_hanging_faces
      end
      hanging_faces_glue = map_parts(hanging_faces_glue...) do hanging_faces_glue...
        hanging_faces_glue
      end
      sDOF_to_dof, sDOF_to_dofs, sDOF_to_coeffs = map_parts(gridap_cell_faces,
                                                            num_regular_faces,
                                                            num_hanging_faces,
                                                            hanging_faces_glue, 
                                                            dmodel.dmodel.models, V.spaces) do gridap_cell_faces,
                                                            num_regular_faces, num_hanging_faces,
                                                            hanging_faces_glue,
                                                            model, V
         
          hanging_faces_to_cell  = Vector{Vector{Int}}(undef, Dc)
          hanging_faces_to_lface = Vector{Vector{Int}}(undef, Dc)

          # Locate for each hanging vertex a cell to which it belongs 
          # and local position within that cell 
          hanging_faces_to_cell[1],    
          hanging_faces_to_lface[1] = _generate_hanging_faces_to_cell_and_lface(num_regular_faces[1], 
                                                                                  num_hanging_faces[1], 
                                                                                  gridap_cell_faces[1])
          
          if (Dc==3)
            hanging_faces_to_cell[2],    
            hanging_faces_to_lface[2] = _generate_hanging_faces_to_cell_and_lface(num_regular_faces[2], 
                                                                                  num_hanging_faces[2], 
                                                                                  gridap_cell_faces[2])
          end 

          # Locate for each hanging facet a cell to which it belongs 
          # and local position within that cell 
          hanging_faces_to_cell[Dc],    
          hanging_faces_to_lface[Dc] = 
                    _generate_hanging_faces_to_cell_and_lface(num_regular_faces[Dc], 
                                                              num_hanging_faces[Dc],
                                                              gridap_cell_faces[Dc])                                                                       

          basis, reffe_args,reffe_kwargs = reffe
          cell_reffe = ReferenceFE(Dc==2 ? QUAD : HEX, basis,reffe_args...;reffe_kwargs...)
          reffe_cell = cell_reffe

          cell_dof_ids  = get_cell_dof_ids(V)
          face_own_dofs = Gridap.ReferenceFEs.get_face_own_dofs(reffe_cell)
          face_dofs     = Gridap.ReferenceFEs.get_face_dofs(reffe_cell)

          hanging_faces_owner_face_dofs = Vector{Vector{Vector{Int}}}(undef, Dc)

          hanging_faces_owner_face_dofs[1] = _generate_hanging_faces_owner_face_dofs(num_hanging_faces[1], 
                                                  face_dofs,
                                                  hanging_faces_glue[1],
                                                  cell_dof_ids)

          if (Dc==3)
            hanging_faces_owner_face_dofs[2] = _generate_hanging_faces_owner_face_dofs(num_hanging_faces[2], 
                                                  face_dofs,
                                                  hanging_faces_glue[2],
                                                  cell_dof_ids)
          end 

          hanging_faces_owner_face_dofs[Dc] = _generate_hanging_faces_owner_face_dofs(num_hanging_faces[Dc], 
                                                  face_dofs,
                                                  hanging_faces_glue[Dc],
                                                  cell_dof_ids)

          sDOF_to_dof    = Int[]
          sDOF_to_dofs   = Vector{Int}[]
          sDOF_to_coeffs = Vector{Float64}[]

          facet_polytope = Dc==2 ? SEGMENT : QUAD
          if (Dc==3)
            edget_polytope = SEGMENT
          end 

          basis, reffe_args,reffe_kwargs = reffe
          face_reffe = ReferenceFE(facet_polytope,basis,reffe_args...;reffe_kwargs...)
          pindex_to_cfvertex_to_fvertex = Gridap.ReferenceFEs.get_vertex_permutations(facet_polytope)

          if (Dc==3)
            edge_reffe = ReferenceFE(edget_polytope,basis,reffe_args...;reffe_kwargs...)
            pindex_to_cevertex_to_evertex = Gridap.ReferenceFEs.get_vertex_permutations(edget_polytope)
          end

          owner_faces_pindex = Vector{Vector{Int}}(undef,Dc-1)
          owner_faces_lids   = Vector{Dict{Int,Tuple{Int,Int,Int}}}(undef, Dc-1)

          lface_to_cvertices = Gridap.ReferenceFEs.get_faces(Dc==2 ? QUAD : HEX, Dc-1, 0)
          owner_faces_pindex[Dc-1], owner_faces_lids[Dc-1]=_compute_owner_faces_pindex_and_lids(Dc-1,Dc,
                                      num_hanging_faces[Dc], 
                                      hanging_faces_glue[Dc],
                                      hanging_faces_to_cell[Dc], 
                                      hanging_faces_to_lface[Dc],
                                      gridap_cell_faces[1],
                                      gridap_cell_faces[Dc],
                                      lface_to_cvertices, 
                                      pindex_to_cfvertex_to_fvertex)

          if (Dc==3)
            ledge_to_cvertices = Gridap.ReferenceFEs.get_faces(HEX, 1, 0)
            pindex_to_cevertex_to_evertex = Gridap.ReferenceFEs.get_vertex_permutations(SEGMENT)
            owner_faces_pindex[1], owner_faces_lids[1]=_compute_owner_faces_pindex_and_lids(1, Dc,
                                      num_hanging_faces[2], 
                                      hanging_faces_glue[2],
                                      hanging_faces_to_cell[2], 
                                      hanging_faces_to_lface[2],
                                      gridap_cell_faces[1],
                                      gridap_cell_faces[2],
                                      ledge_to_cvertices, 
                                      pindex_to_cevertex_to_evertex)
          end
        
  
          node_permutations = Vector{Vector{Vector{Int}}}(undef, Dc-1)
          nodes, _ = Gridap.ReferenceFEs.compute_nodes(facet_polytope,[reffe_args[2] for i=1:Dc-1])
          node_permutations[Dc-1]=Gridap.ReferenceFEs._compute_node_permutations(facet_polytope,nodes)
          if (Dc==3)
            nodes, _ = Gridap.ReferenceFEs.compute_nodes(edget_polytope,[reffe_args[2] for i=1:Dc-2])
            node_permutations[1]=Gridap.ReferenceFEs._compute_node_permutations(edget_polytope,nodes)
          end 
          
          subface_own_dofs = Vector{Vector{Vector{Int}}}(undef, Dc-1)
          subface_own_dofs[Dc-1] = Gridap.ReferenceFEs.get_face_own_dofs(face_reffe)
          if (Dc==3)
            subface_own_dofs[1] = Gridap.ReferenceFEs.get_face_own_dofs(edge_reffe)
          end 
          _generate_constraints!(0,
                                Dc,
                                [gridap_cell_faces[i] for i=1:Dc],
                                num_hanging_faces[1],
                                hanging_faces_to_cell[1],
                                hanging_faces_to_lface[1],
                                hanging_faces_owner_face_dofs[1],
                                hanging_faces_glue[1],
                                face_subface_ldof_to_cell_ldof,
                                face_dofs,
                                face_own_dofs,
                                subface_own_dofs,
                                cell_dof_ids,
                                node_permutations,
                                owner_faces_pindex,
                                owner_faces_lids,
                                ref_constraints,
                                sDOF_to_dof,
                                sDOF_to_dofs,
                                sDOF_to_coeffs)
          
          _generate_constraints!(Dc-1,
                                Dc,
                                [gridap_cell_faces[i] for i=1:Dc],
                                num_hanging_faces[Dc],
                                hanging_faces_to_cell[Dc],
                                hanging_faces_to_lface[Dc],
                                hanging_faces_owner_face_dofs[Dc],
                                hanging_faces_glue[Dc],
                                face_subface_ldof_to_cell_ldof,
                                face_dofs,
                                face_own_dofs,
                                subface_own_dofs,
                                cell_dof_ids,
                                node_permutations,
                                owner_faces_pindex,
                                owner_faces_lids,
                                ref_constraints,
                                sDOF_to_dof,
                                sDOF_to_dofs,
                                sDOF_to_coeffs)
          sDOF_to_dof, Gridap.Arrays.Table(sDOF_to_dofs), Gridap.Arrays.Table(sDOF_to_coeffs)
        end 
    end

    rr = Gridap.Adaptivity.RedRefinementRule(cell_polytope)
    
    face_subface_ldof_to_cell_ldof = Vector{Vector{Vector{Vector{Int32}}}}(undef,Dc-1)
    face_subface_ldof_to_cell_ldof[Dc-1] = 
        Gridap.Adaptivity.get_face_subface_ldof_to_cell_ldof(rr,Tuple(order for _=1:Dc),Dc-1)
    if (Dc==3)
      face_subface_ldof_to_cell_ldof[1] =
        Gridap.Adaptivity.get_face_subface_ldof_to_cell_ldof(rr,Tuple(order for _=1:Dc),1)
    end 

    sDOF_to_dof, sDOF_to_dofs,sDOF_to_coeffs=
        generate_constraints(dmodel, V, reffe, 
                            non_conforming_glue, ref_constraints, face_subface_ldof_to_cell_ldof)
    println(sDOF_to_dof)
    println(sDOF_to_dofs)
    println(sDOF_to_coeffs)

    # Define manufactured functions
    u(x) = x[1]+x[2]^order
    f(x) = -Δ(u)(x)

    map_parts(dmodel.dmodel.models,V.spaces,U.spaces,sDOF_to_dof,sDOF_to_dofs,sDOF_to_coeffs) do model,V,U,sDOF_to_dof,sDOF_to_dofs,sDOF_to_coeffs
      
      println(get_cell_dof_ids(V))
      fl=get_face_labeling(model)
      t=GridapDistributed.get_grid_topology(model)
      println(Gridap.Geometry.get_faces(t,2,0))

      for i=1:length(fl.d_to_dface_to_entity)-1
        println(fl.d_to_dface_to_entity[i])
      end  

      Vc = FESpaceWithLinearConstraints(
        sDOF_to_dof,
        sDOF_to_dofs,
        sDOF_to_coeffs,
        V)
      Uc = TrialFESpace(Vc,u)

      # Define integration mesh and quadrature
      degree = 2*order+1
      Ω = Triangulation(model)
      dΩ = Measure(Ω,degree)

      a(u,v) = ∫( ∇(v)⊙∇(u) )*dΩ
      b(v) = ∫(v*f)*dΩ

      # op = AffineFEOperator(a,b,U,V0)
      op = AffineFEOperator(a,b,Uc,Vc)
      uh = solve(op)

      # Define exact solution and error
      e = u - uh

      # Compute errors
      el2 = sqrt(sum( ∫( e*e )*dΩ ))
      eh1 = sqrt(sum( ∫( e*e + ∇(e)⋅∇(e) )*dΩ ))

      tol=1e-8
      @assert el2 < tol
      @assert eh1 < tol
    end 
  end 

  test(Val{2},1,1)
  test(Val{2},1,2)
  test(Val{2},1,3)

  test(Val{2},2,1)
  test(Val{2},2,2)
  test(Val{2},2,3)

  test(Val{2},3,1)
  test(Val{2},3,2)
  test(Val{2},3,3)

  test(Val{2},4,1)
  test(Val{2},4,2)
  test(Val{2},4,3)

  test(Val{3},1,1)
  test(Val{3},1,2)

#end