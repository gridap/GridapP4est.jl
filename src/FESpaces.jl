
function _build_constraint_coefficients_matrix_in_ref_space(Dc, reffe::Tuple{<:Lagrangian,Any,Any})
    cell_polytope = Dc == 2 ? QUAD : HEX
    basis, reffe_args, reffe_kwargs = reffe
    cell_reffe = ReferenceFE(cell_polytope, basis, reffe_args...; reffe_kwargs...)
    
    # TO-DO: How can modelH be created such that it is tailored to cell_polytope?
    modelH= _generate_unit_hypercube_model(Dc)
    modelh=refine(modelH,2)
  
    VH=TestFESpace(modelH,cell_reffe)
    Vh=TestFESpace(modelh,cell_reffe)

    uH=get_fe_basis(VH)
    uHh=change_domain(uH,get_triangulation(modelh),ReferenceDomain())
    σRTh=Gridap.FESpaces.get_fe_dof_basis(Vh)
    ref_constraints_contribs=σRTh(uHh) 

    ref_constraints = Matrix{Float64}(undef,num_free_dofs(Vh),num_free_dofs(VH))
    cell_dof_ids = get_cell_dof_ids(Vh)
    cache_cell_dof_ids = array_cache(cell_dof_ids)
    cache_ref_constraints_contribs = array_cache(ref_constraints_contribs) 
    for cell=1:length(cell_dof_ids)
       current_cell_dof_ids=getindex!(cache_cell_dof_ids,cell_dof_ids,cell)
       current_ref_constraints_contribs=getindex!(cache_ref_constraints_contribs,ref_constraints_contribs,cell)
       ref_constraints[current_cell_dof_ids,:]=current_ref_constraints_contribs
    end
    ref_constraints
end

function _fill_face_subface_ldof_to_cell_ldof!(face_subface_ldof_to_cell_ldof,
                                               num_faces,
                                               coarse_faces_to_child_ids,
                                               face_dofs,
                                               cells_dof_ids,
                                               first_face)
    for coarse_face_id=1:num_faces 
        for (subface,child_id) in enumerate(coarse_faces_to_child_ids[coarse_face_id,:])
            @debug "coarse_face_id: $(coarse_face_id), subface: $(subface), child_id: $(child_id)"
            cell_dof_ids=cells_dof_ids[child_id]
            for (i,dof) in enumerate(cell_dof_ids[face_dofs[first_face+coarse_face_id]])
                @debug "i: $(i), dof: $(dof)"
                face_subface_ldof_to_cell_ldof[coarse_face_id][subface][i]=dof 
            end
        end 
    end
end 


function _generate_face_subface_ldof_to_cell_ldof(Df,Dc,reffe::Tuple{<:Lagrangian,Any,Any})
    cell_polytope = (Dc == 2) ? QUAD : HEX
    coarse_faces_to_child_ids = (Dc == 2) ? _coarse_faces_to_child_ids_2D : 
                                            _coarse_faces_to_child_ids_3D

    basis, reffe_args, reffe_kwargs = reffe
    cell_reffe = ReferenceFE(cell_polytope, basis, reffe_args...; reffe_kwargs...)

    # TO-DO: How can modelH be created such that it is tailored to cell_polytope?
    modelH= _generate_unit_hypercube_model(Dc)
    modelh=refine(modelH,2)
    Vh=TestFESpace(modelh,reffe)
    cells_dof_ids=get_cell_dof_ids(Vh)

    first_face = get_offset(get_polytope(cell_reffe),Df)
    face_dofs=get_face_dofs(cell_reffe)
    num_dofs_x_face = length(face_dofs[first_face+1])
    
    if (Df==Dc-1) # Facets    
        num_faces = 2*Dc
        num_subfaces = 2^(Dc-1)
        face_subface_ldof_to_cell_ldof=_allocate_face_subface_ldof_to_cell_ldof(num_faces,num_subfaces,num_dofs_x_face)
        _fill_face_subface_ldof_to_cell_ldof!(face_subface_ldof_to_cell_ldof,
                                               num_faces,
                                               coarse_faces_to_child_ids,
                                               face_dofs,
                                               cells_dof_ids,
                                               first_face)
        face_subface_ldof_to_cell_ldof
    else
        @assert Df==1 && Dc==3 # Edges in 3D
        num_edges=12
        num_subedges=2
        edge_subedge_ldof_to_cell_ldof=_allocate_face_subface_ldof_to_cell_ldof(num_edges,num_subedges,num_dofs_x_face)
        _fill_face_subface_ldof_to_cell_ldof!(edge_subedge_ldof_to_cell_ldof,
                                             num_edges,
                                              _coarse_edges_to_child_ids_3D,
                                              face_dofs,
                                              cells_dof_ids,
                                              first_face)
        edge_subedge_ldof_to_cell_ldof
    end 
end 



const _coarse_faces_to_child_ids_2D=[1 2; 3 4; 1 3; 2 4]
const _coarse_faces_to_child_ids_3D=[1 2 3 4; 5 6 7 8; 1 2 5 6; 3 4 7 8; 1 3 5 7; 2 4 6 8; ]
const _coarse_edges_to_child_ids_3D=[1 2; 3 4; 5 6; 7 8; 1 3; 2 4; 5 7; 6 8; 1 5; 2 6; 3 7; 4 8; ]


function _generate_unit_hypercube_model(Dc)
    @assert Dc==2 || Dc==3
    if (Dc==2)
      modelH=CartesianDiscreteModel((0,1,0,1),(1,1))
    else
      modelH=CartesianDiscreteModel((0,1,0,1,0,1),(1,1,1))   
    end
    modelH
end 

function _allocate_face_subface_ldof_to_cell_ldof(num_faces, num_subfaces, num_dofs_x_face)
    face_subface_ldof_to_cell_ldof = Vector{Vector{Vector{Int}}}(undef,num_faces)
    for face=1:num_faces
        face_subface_ldof_to_cell_ldof[face]=Vector{Vector{Int}}(undef,num_subfaces)
        for subface=1:num_subfaces
            face_subface_ldof_to_cell_ldof[face][subface]=Vector{Int}(undef,num_dofs_x_face)
        end
    end 
    face_subface_ldof_to_cell_ldof
end 

function _generate_face_subface_ldof_to_cell_ldof(Df,Dc,reffe::Tuple{<:RaviartThomas,Any,Any})
    cell_polytope = (Dc == 2) ? QUAD : HEX
    coarse_faces_to_child_ids = (Dc == 2) ? _coarse_faces_to_child_ids_2D : 
                                            _coarse_faces_to_child_ids_3D
    
    if (Df==Dc-1) # Facets    
        basis, reffe_args, reffe_kwargs = reffe
        cell_reffe = ReferenceFE(cell_polytope, basis, reffe_args...; reffe_kwargs...)

        # TO-DO: How can modelH be created such that it is tailored to cell_polytope?
        modelH= _generate_unit_hypercube_model(Dc)
        modelh=refine(modelH,2)
        RTh=TestFESpace(modelh,reffe)
        num_faces = 2*Dc
        num_subfaces = 2^(Dc-1)
        first_face = get_offset(get_polytope(cell_reffe),Df)
        face_own_dofs=get_face_own_dofs(cell_reffe)
        num_dofs_x_face = length(face_own_dofs[first_face+1])
        face_subface_ldof_to_cell_ldof=_allocate_face_subface_ldof_to_cell_ldof(num_faces,num_subfaces,num_dofs_x_face)
        RTh_cell_dof_ids=get_cell_dof_ids(RTh)
        _fill_face_subface_ldof_to_cell_ldof!(face_subface_ldof_to_cell_ldof,
                                        num_faces,
                                        coarse_faces_to_child_ids,
                                        face_own_dofs,
                                        RTh_cell_dof_ids,
                                        first_face)
    else
        @assert Df==1 && Dc==3 # Edges in 3D 
        num_edges=12
        num_subedges=2 
        num_dofs_x_edge=0
        face_subface_ldof_to_cell_ldof=_allocate_face_subface_ldof_to_cell_ldof(num_edges,
                                                                                num_subedges,
                                                                                num_dofs_x_edge)        
    end 
    face_subface_ldof_to_cell_ldof
end 

function _build_constraint_coefficients_matrix_in_ref_space(Dc, reffe::Tuple{<:RaviartThomas,Any,Any})
    cell_polytope = Dc == 2 ? QUAD : HEX
    basis, reffe_args, reffe_kwargs = reffe
    cell_reffe = ReferenceFE(cell_polytope, basis, reffe_args...; reffe_kwargs...)
    
    # TO-DO: How can modelH be created such that it is tailored to cell_polytope?
    modelH= _generate_unit_hypercube_model(Dc)
    modelh=refine(modelH,2)
  
    VH=TestFESpace(modelH,cell_reffe)
    Vh=TestFESpace(modelh,cell_reffe)

    uH=get_fe_basis(VH)
    uHh=change_domain(uH,get_triangulation(modelh),ReferenceDomain())
    σRTh=Gridap.FESpaces.get_fe_dof_basis(Vh)
    ref_constraints_contribs=σRTh(uHh) 

    ref_constraints = Matrix{Float64}(undef,num_free_dofs(Vh),num_free_dofs(VH))
    cell_dof_ids = get_cell_dof_ids(Vh)
    cache_cell_dof_ids = array_cache(cell_dof_ids)
    cache_ref_constraints_contribs = array_cache(ref_constraints_contribs) 
    for cell=1:length(cell_dof_ids)
       current_cell_dof_ids=getindex!(cache_cell_dof_ids,cell_dof_ids,cell)
       current_ref_constraints_contribs=getindex!(cache_ref_constraints_contribs,ref_constraints_contribs,cell)
       ref_constraints[current_cell_dof_ids,:]=current_ref_constraints_contribs
    end
    # We change the sign of the constraints coefficients here so that 
    # the unit normal to hanging faces matches the unit normal 
    # of the owner face. This way, we properly glue the global DoFs at 
    # both sides of the interface of cells at different refinement level.
    # We could have done this instead by adjusting the sign_flips for hanging
    # faces in the FESpace constructor, although we decide it to do it here 
    # because it is way simpler.
    @. ref_constraints = -ref_constraints
    ref_constraints
end


# To-think: might this info go to the glue? 
# If it is required in different scenarios, I would say it may make sense
function _generate_hanging_faces_to_cell_and_lface(num_regular_faces,
    num_hanging_faces,
    gridap_cell_faces)
    # Locate for each hanging face the owner cell among the set of cells 
    # to which it belongs and local position within that cell 
    # By convention, the owner cell is the cell with minimum identifer among 
    # all in the set. This convention is used here, and in other parts of the code.
    # Breaking this convention here may affect the consistency of other parts of the code. 
    hanging_faces_to_cell = Vector{Int}(undef, num_hanging_faces)
    hanging_faces_to_lface = Vector{Int}(undef, num_hanging_faces)
    hanging_faces_to_cell .= -1 
    for cell = 1:length(gridap_cell_faces)
        s = gridap_cell_faces.ptrs[cell]
        e = gridap_cell_faces.ptrs[cell+1]
        l = e - s
        for j = 1:l
            fid = gridap_cell_faces.data[s+j-1]
            if fid > num_regular_faces
                fid_hanging = fid - num_regular_faces
                if (hanging_faces_to_cell[fid_hanging]==-1)
                  hanging_faces_to_cell[fid_hanging] = cell
                  hanging_faces_to_lface[fid_hanging] = j
                end
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
    ptrs = Vector{Int}(undef, num_hanging_faces + 1)
    ptrs[1] = 1
    for fid_hanging = 1:num_hanging_faces
        glue = hanging_faces_glue[fid_hanging]
        if (glue[1]!=-1)
          ocell_lface = glue[2]
          ptrs[fid_hanging+1] = ptrs[fid_hanging] + length(face_dofs[ocell_lface])
        else
          ptrs[fid_hanging+1] = ptrs[fid_hanging] + 1 
        end 
    end
    data_owner_face_dofs = Vector{Int}(undef, ptrs[num_hanging_faces+1] - 1)
    for fid_hanging = 1:num_hanging_faces
        glue = hanging_faces_glue[fid_hanging]
        ocell = glue[1]
        if (ocell!=-1)
          ocell_lface = glue[2]
          s = ptrs[fid_hanging]
          e = ptrs[fid_hanging+1] - 1
          current_cell_dof_ids = getindex!(cache, cell_dof_ids, ocell)
          for (j, ldof) in enumerate(face_dofs[ocell_lface])
            data_owner_face_dofs[s+j-1] = current_cell_dof_ids[ldof]
          end
        else 
          s = ptrs[fid_hanging]
          data_owner_face_dofs[s] = 1
        end
    end
    Gridap.Arrays.Table(data_owner_face_dofs, ptrs)
end

function face_dim(::Type{Val{Dc}}, face_lid) where {Dc}
    num_vertices = GridapP4est.num_cell_vertices(Val{Dc})
    num_edges = GridapP4est.num_cell_edges(Val{Dc})
    num_faces = GridapP4est.num_cell_faces(Val{Dc})
    if (face_lid <= num_vertices)
        return 0
    elseif (face_lid <= num_vertices + num_edges)
        return 1
    elseif (face_lid <= num_vertices + num_edges + num_faces)
        return Dc - 1
    end
end

function face_lid_within_dim(::Type{Val{Dc}}, face_lid) where {Dc}
    num_vertices = GridapP4est.num_cell_vertices(Val{Dc})
    num_edges = GridapP4est.num_cell_edges(Val{Dc})
    num_faces = GridapP4est.num_cell_faces(Val{Dc})
    if (face_lid <= num_vertices)
        return face_lid
    elseif (face_lid <= num_vertices + num_edges)
        return face_lid - num_vertices
    elseif (face_lid <= num_vertices + num_edges + num_faces)
        return face_lid - num_vertices - num_edges
    end
end

function _restrict_face_dofs_to_face_dim(cell_reffe,Df)
    polytope = get_polytope(cell_reffe)
    first_face = Gridap.ReferenceFEs.get_offset(polytope,Df)+1
    face_own_dofs=get_face_own_dofs(cell_reffe)
    first_face_faces = Gridap.ReferenceFEs.get_faces(polytope)[first_face]
    face_dofs_to_face_dim = face_own_dofs[first_face_faces]
    touched=Dict{Int,Int}()
    dofs=Int64[]
    current=1
    for (i,face_dofs) in enumerate(face_dofs_to_face_dim)
        for (j,dof) in enumerate(face_dofs)
            push!(dofs,dof)
        end
    end
    sort!(dofs)
    touched = Dict( dofs[i]=>i for i=1:length(dofs))
    for (i,face_dofs) in enumerate(face_dofs_to_face_dim)
        for (j,dof) in enumerate(face_dofs)
            face_dofs[j]=touched[dof]
        end
    end
    face_dofs_to_face_dim
end 
function _num_face_own_dofs(cell_reffe,Df)
  polytope = get_polytope(cell_reffe)
  first_face = Gridap.ReferenceFEs.get_offset(polytope,Df)+1
  face_own_dofs=get_face_own_dofs(cell_reffe)
  length(face_own_dofs[first_face])
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

    @assert Dc == 2 || Dc == 3
    @assert 0 ≤ Df < Dc

    num_vertices = GridapP4est.num_cell_vertices(Val{Dc})
    num_edges = GridapP4est.num_cell_edges(Val{Dc})
    num_faces = GridapP4est.num_cell_faces(Val{Dc})

    offset = 0
    if (Df ≥ 1)
        offset += num_vertices
        if (Df == 2)
            offset += num_edges
        end
    end

    cache_dof_ids = array_cache(cell_dof_ids)
    first_cache_cell_faces = array_cache(first(cell_faces))
    cache_cell_faces = Vector{typeof(first_cache_cell_faces)}(undef, length(cell_faces))
    for i = 1:length(cell_faces)
        cache_cell_faces[i] = array_cache(cell_faces[i])
    end

    for fid_hanging = 1:num_hanging_faces
        cell = hanging_faces_to_cell[fid_hanging]
        current_cell_dof_ids = getindex!(cache_dof_ids, cell_dof_ids, cell)
        lface = hanging_faces_to_lface[fid_hanging]
        ocell, ocell_lface, subface = hanging_faces_glue[fid_hanging]
        if (ocell!=-1)
            ocell_lface_within_dim = face_lid_within_dim(Val{Dc}, ocell_lface)
            oface_dim = face_dim(Val{Dc}, ocell_lface)

            if (Df == 0) # Am I a vertex?
                hanging_lvertex_within_first_subface = 2^oface_dim
                cur_subface_own_dofs = subface_own_dofs[oface_dim][hanging_lvertex_within_first_subface]
            elseif (Df == 1 && Dc == 3) # Am I an edge?
                if (subface < 0) # Edge hanging in the interior of a face 
                    @assert subface == -1 || subface == -2 || subface == -3 || subface == -4
                    @assert oface_dim == Dc - 1
                    abs_subface = abs(subface)
                    if (abs_subface == 1)
                        subface = 1
                        edge = 4 + 4 # num_vertices+edge_id
                    elseif (abs_subface == 2)
                        subface = 3
                        edge = 4 + 4 # num_vertices+edge_id 
                    elseif (abs_subface == 3)
                        subface = 3
                        edge = 4 + 1 # num_vertices+edge_id 
                    elseif (abs_subface == 4)
                        subface = 4
                        edge = 4 + 1 # num_vertices+edge_id 
                    end
                    cur_subface_own_dofs = subface_own_dofs[oface_dim][edge]
                else
                    @assert subface == 1 || subface == 2
                    @assert oface_dim == 1
                    hanging_lvertex_within_first_subface = 2^oface_dim
                    cur_subface_own_dofs = subface_own_dofs[oface_dim][end]
                end
            elseif (Df == Dc - 1) # Am I a face?
                @assert oface_dim == Dc - 1
                cur_subface_own_dofs = subface_own_dofs[oface_dim][end]
            end


            oface = getindex!(cache_cell_faces[oface_dim+1], cell_faces[oface_dim+1], ocell)[ocell_lface_within_dim]
            oface_lid, _ = owner_faces_lids[oface_dim][oface]
            pindex = owner_faces_pindex[oface_dim][oface_lid]
            @debug "Df=$(Df) Dc=$(Dc) fid_hanging=$(fid_hanging) cell=$(cell) oface=$(oface) oface_lid=$(oface_lid) pindex=$(pindex) ocell=$(ocell) ocell_lface=$(ocell_lface) subface=$(subface) oface_dim=$(oface_dim) cur_subface_own_dofs=$(cur_subface_own_dofs) face_own_dofs=$(face_own_dofs[offset+lface])"
            for ((ldof, dof), ldof_subface) in zip(enumerate(face_own_dofs[offset+lface]), cur_subface_own_dofs)
                push!(sDOF_to_dof, current_cell_dof_ids[dof])
                push!(sDOF_to_dofs, hanging_faces_owner_face_dofs[fid_hanging])
                coeffs = Vector{Float64}(undef, length(hanging_faces_owner_face_dofs[fid_hanging]))
                # Go over dofs of ocell_lface
                for (ifdof, icdof) in enumerate(face_dofs[ocell_lface])
                    pifdof = node_permutations[oface_dim][pindex][ifdof]
                    ldof_coarse = face_dofs[ocell_lface][pifdof]
                    coeffs[ifdof] =
                        ref_constraints[face_subface_ldof_to_cell_ldof[oface_dim][ocell_lface_within_dim][subface][ldof_subface], ldof_coarse]
                end
                push!(sDOF_to_coeffs, coeffs)
            end
        else 
            for (ldof, dof) in enumerate(face_own_dofs[offset+lface])
                push!(sDOF_to_dof, current_cell_dof_ids[dof])
                push!(sDOF_to_dofs, hanging_faces_owner_face_dofs[fid_hanging])
                push!(sDOF_to_coeffs, [1.0])
            end
        end 
    end
    @debug "sDOF_to_dof [$(Df)]= $(sDOF_to_dof)"
    @debug "sDOF_to_dofs [$(Df)]= $(sDOF_to_dofs)"
    @debug "sDOF_to_coeffs [$(Df)]= $(sDOF_to_coeffs)"
end

function _compute_owner_faces_lids(Df,Dc,num_hanging_faces,hanging_faces_glue,cell_faces)
    num_owner_faces = 0
    owner_faces_lids = Dict{Int,Tuple{Int,Int,Int}}()
    for fid_hanging = 1:num_hanging_faces
        ocell, ocell_lface, _ = hanging_faces_glue[fid_hanging]
        if (ocell!=-1)
            ocell_dim = face_dim(Val{Dc}, ocell_lface)
            if (ocell_dim == Df)
                ocell_lface_within_dim = face_lid_within_dim(Val{Dc}, ocell_lface)
                owner_face = cell_faces[ocell][ocell_lface_within_dim]
                if !(haskey(owner_faces_lids, owner_face))
                    num_owner_faces += 1
                    owner_faces_lids[owner_face] = (num_owner_faces, ocell, ocell_lface)
                end
            end
        end
    end
    owner_faces_lids
end 


# count how many different owner faces
# for each owner face 
#    track the global IDs of its face vertices from the perspective of the subfaces
# for each owner face 
#    compute permutation id
function _compute_owner_faces_pindex_and_lids(Dc,
                                              num_hanging_faces,
                                              hanging_faces_glue,
                                              hanging_faces_to_cell,
                                              hanging_faces_to_lface,
                                              cell_vertices,
                                              cell_faces,
                                              lface_to_cvertices,
                                              pindex_to_cfvertex_to_fvertex)

    owner_faces_lids =_compute_owner_faces_lids(Dc-1,Dc,
                      num_hanging_faces,hanging_faces_glue,cell_faces)                                       
    
    @debug "owner_faces_lids [Df=$(Dc-1) Dc=$(Dc)]: $(owner_faces_lids)"

    num_owner_faces   = length(keys(owner_faces_lids))
    num_face_vertices = length(first(lface_to_cvertices))
    owner_face_vertex_ids = Vector{Int}(undef, num_face_vertices * num_owner_faces)
    owner_face_vertex_ids .= -1

    for fid_hanging = 1:num_hanging_faces
        ocell, ocell_lface, subface = hanging_faces_glue[fid_hanging]
        @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))]: fid_hanging=$(fid_hanging) ocell=$(ocell) ocell_lface=$(ocell_lface) subface=$(subface)"

        if (ocell!=-1)
            oface_dim = face_dim(Val{Dc}, ocell_lface)
            if (oface_dim == Dc-1)
                cell = hanging_faces_to_cell[fid_hanging]
                lface = hanging_faces_to_lface[fid_hanging]
                cvertex = lface_to_cvertices[lface][subface]
                vertex = cell_vertices[cell][cvertex]
                ocell_lface_within_dim = face_lid_within_dim(Val{Dc}, ocell_lface)
                owner_face = cell_faces[ocell][ocell_lface_within_dim]
                owner_face_lid, _ = owner_faces_lids[owner_face]
                @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))]: cell=$(cell) lface=$(lface) cvertex=$(cvertex) vertex=$(vertex) owner_face=$(owner_face) owner_face_lid=$(owner_face_lid)"
                @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))]: owner_face_vertex_ids[$((owner_face_lid-1)*num_face_vertices+subface)] = $(vertex)"
                owner_face_vertex_ids[(owner_face_lid-1)*num_face_vertices+subface] = vertex
            end
        end
    end
    @debug "owner_face_vertex_ids [Dc=$(Dc)]: $(owner_face_vertex_ids)"

    owner_faces_pindex = Vector{Int}(undef, num_owner_faces)
    for owner_face in keys(owner_faces_lids)
        (owner_face_lid, ocell, ocell_lface) = owner_faces_lids[owner_face]
        ocell_lface_within_dim = face_lid_within_dim(Val{Dc}, ocell_lface)
        # Compute permutation id by comparing 
        #  1. cell_vertices[ocell][ocell_lface]
        #  2. owner_face_vertex_ids 
        pindexfound = false
        cfvertex_to_cvertex = lface_to_cvertices[ocell_lface_within_dim]
        for (pindex, cfvertex_to_fvertex) in enumerate(pindex_to_cfvertex_to_fvertex)
            found = true
            for (cfvertex, fvertex) in enumerate(cfvertex_to_fvertex)
                vertex1 = owner_face_vertex_ids[(owner_face_lid-1)*num_face_vertices+fvertex]
                cvertex = cfvertex_to_cvertex[cfvertex]
                vertex2 = cell_vertices[ocell][cvertex]
                @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))]: cell_vertices[$(ocell)][$(cvertex)]=$(cell_vertices[ocell][cvertex]) owner_face_vertex_ids[$((owner_face_lid-1)*num_face_vertices+fvertex)]=$(owner_face_vertex_ids[(owner_face_lid-1)*num_face_vertices+fvertex])"
                # -1 can only happen in the interface of two 
                # ghost cells at different refinement levels
                if (vertex1 != vertex2) && (vertex1 != -1) 
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
    @debug "owner_faces_pindex: $(owner_faces_pindex)"

    owner_faces_pindex, owner_faces_lids
end


function _compute_owner_edges_pindex_and_lids(
    num_hanging_edges,
    hanging_edges_glue,
    hanging_edges_to_cell,
    hanging_edges_to_ledge,
    cell_vertices,
    cell_edges)
    Dc=3
    owner_edges_lids =_compute_owner_faces_lids(1,
                                                Dc,
                                                num_hanging_edges,
                                                hanging_edges_glue,
                                                cell_edges)

    num_owner_edges = length(keys(owner_edges_lids))
    owner_edges_pindex = Vector{Int}(undef, num_owner_edges)

    ledge_to_cvertices = Gridap.ReferenceFEs.get_faces(HEX, 1, 0)

    # Go over hanging edges 
    # Find the owner hanging edge
    for fid_hanging = 1:num_hanging_edges
        ocell, ocell_ledge, subedge = hanging_edges_glue[fid_hanging]
        ocell_dim = face_dim(Val{Dc}, ocell_ledge)
        if (ocell!=-1 && ocell_dim==1)
          ocell_ledge_within_dim = face_lid_within_dim(Val{Dc}, ocell_ledge)
          cell = hanging_edges_to_cell[fid_hanging]
          ledge = hanging_edges_to_ledge[fid_hanging]
          gvertex1 = cell_vertices[cell][ledge_to_cvertices[ledge][subedge]]
          gvertex2 = cell_vertices[ocell][ledge_to_cvertices[ocell_ledge_within_dim][subedge]]
          @debug "fid_hanging=$(fid_hanging) cell=$(cell) ledge=$(ledge) ocell=$(ocell) ocell_ledge=$(ocell_ledge) subedge=$(subedge) gvertex1=$(gvertex1) gvertex2=$(gvertex2)"
          pindex = gvertex1==gvertex2 ? 1 : 2
          owner_edge=cell_edges[ocell][ocell_ledge_within_dim]
          owner_edge_lid, _ = owner_edges_lids[owner_edge]
          owner_edges_pindex[owner_edge_lid]=pindex
        end
    end
    owner_edges_pindex, owner_edges_lids
end

using Gridap.ReferenceFEs

function get_nodes_permutations(reffe::GenericLagrangianRefFE{GradConformity})
  p = get_polytope(reffe)
  face_nodes = get_face_nodes(reffe)
  dofs = get_dof_basis(reffe)
  interior_nodes = dofs.nodes[face_nodes[end]]
  Gridap.ReferenceFEs._compute_node_permutations(p,interior_nodes)
end

function get_nodes_permutations(reffe::ReferenceFE,d::Integer)
  p = get_polytope(reffe)
  range = get_dimrange(p,d)
  get_face_nodes_permutations(reffe)[range]
end

function get_face_nodes_permutations(reffe::GenericLagrangianRefFE{GradConformity})
    nodes_permutations = get_nodes_permutations(reffe)
    reffaces = reffe.reffe.metadata
    _reffaces = vcat(reffaces...)
    face_nodes_permutations = map(get_nodes_permutations,_reffaces)
    push!(face_nodes_permutations,nodes_permutations)
    face_nodes_permutations
end

function get_face_dofs_permutations(reffe::LagrangianRefFE)
    dofs = get_dof_basis(reffe)
    face_nodes_permutations = get_face_nodes_permutations(reffe)
    face_nodes = get_face_nodes(reffe)
    face_dofs = get_face_dofs(reffe)
    face_dofs_permutations = Gridap.ReferenceFEs._generate_face_own_dofs_permutations(
      face_nodes_permutations, dofs.node_and_comp_to_dof, face_nodes, face_dofs)
    face_dofs_permutations
end

function get_face_dofs_permutations(reffe::ReferenceFE,d::Integer)
    p = get_polytope(reffe)
    range = get_dimrange(p,d)
    get_face_dofs_permutations(reffe)[range]
end

function get_face_dofs_permutations(
         reffe::Gridap.ReferenceFEs.GenericRefFE{Gridap.ReferenceFEs.RaviartThomas, Dc},Df::Integer) where Dc
    first_face = get_offset(get_polytope(reffe),Df)
    order = length(get_face_dofs(reffe)[first_face])-1
    nfaces=num_faces(reffe,Df)
    if (Df==Dc-1)
       facet_polytope = Dc == 2 ? SEGMENT : QUAD
       nodes, _ = Gridap.ReferenceFEs.compute_nodes(facet_polytope, [order for i = 1:Dc-1])
       Fill(Gridap.ReferenceFEs._compute_node_permutations(facet_polytope, nodes),nfaces)
    elseif (Dc == 3 && Df==1)
       nodes, _ = Gridap.ReferenceFEs.compute_nodes(SEGMENT, [order for i = 1:Dc-2])
       Fill(Gridap.ReferenceFEs._compute_node_permutations(SEGMENT, nodes),nfaces)
    end
end 

function generate_constraints(dmodel::OctreeDistributedDiscreteModel{Dc},
    spaces_wo_constraints,
    reffe,
    ref_constraints,
    face_subface_ldof_to_cell_ldof) where {Dc}

    non_conforming_glue = dmodel.non_conforming_glue
    dmodel = dmodel.dmodel

    gridap_cell_faces = map(local_views(dmodel)) do model
        topo = Gridap.Geometry.get_grid_topology(model)
        Tuple(Gridap.Geometry.get_faces(topo, Dc, d) for d = 0:Dc-1)
    end
    num_regular_faces = map(non_conforming_glue) do ncglue
        @debug "num_regular_faces=$(Tuple(ncglue.num_regular_faces[d] for d = 1:Dc))"
        Tuple(ncglue.num_regular_faces[d] for d = 1:Dc)
    end
    num_hanging_faces = map(non_conforming_glue) do ncglue
        @debug "num_hanging_faces=$(Tuple(ncglue.num_hanging_faces[d] for d = 1:Dc))"
        Tuple(ncglue.num_hanging_faces[d] for d = 1:Dc)
    end
    hanging_faces_glue = map(non_conforming_glue) do ncglue
        Tuple(ncglue.hanging_faces_glue[d] for d = 1:Dc)
    end
    sDOF_to_dof, sDOF_to_dofs, sDOF_to_coeffs = map(gridap_cell_faces,
        num_regular_faces,
        num_hanging_faces,
        hanging_faces_glue,
        dmodel.models,
        spaces_wo_constraints) do gridap_cell_faces,
                                    num_regular_faces,
                                    num_hanging_faces,
                                    hanging_faces_glue,
                                    model,
                                    V

        hanging_faces_to_cell = Vector{Vector{Int}}(undef, Dc)
        hanging_faces_to_lface = Vector{Vector{Int}}(undef, Dc)

        # Locate for each hanging vertex a cell to which it belongs 
        # and local position within that cell 
        hanging_faces_to_cell[1],
        hanging_faces_to_lface[1] = _generate_hanging_faces_to_cell_and_lface(num_regular_faces[1],
            num_hanging_faces[1],
            gridap_cell_faces[1])

        if (Dc == 3)
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

        basis, reffe_args, reffe_kwargs = reffe
        cell_reffe = ReferenceFE(Dc == 2 ? QUAD : HEX, basis, reffe_args...; reffe_kwargs...)

        cell_dof_ids = get_cell_dof_ids(V)
        face_own_dofs = Gridap.ReferenceFEs.get_face_own_dofs(cell_reffe)
        face_dofs = Gridap.ReferenceFEs.get_face_dofs(cell_reffe)

        hanging_faces_owner_face_dofs = Vector{Vector{Vector{Int}}}(undef, Dc)

        hanging_faces_owner_face_dofs[1] = _generate_hanging_faces_owner_face_dofs(num_hanging_faces[1],
            face_dofs,
            hanging_faces_glue[1],
            cell_dof_ids)

        if (Dc == 3)
            hanging_faces_owner_face_dofs[2] = _generate_hanging_faces_owner_face_dofs(num_hanging_faces[2],
                face_dofs,
                hanging_faces_glue[2],
                cell_dof_ids)
        end

        hanging_faces_owner_face_dofs[Dc] = _generate_hanging_faces_owner_face_dofs(num_hanging_faces[Dc],
            face_dofs,
            hanging_faces_glue[Dc],
            cell_dof_ids)

        sDOF_to_dof = Int[]
        sDOF_to_dofs = Vector{Int}[]
        sDOF_to_coeffs = Vector{Float64}[]

        facet_polytope = Dc == 2 ? SEGMENT : QUAD
        if (Dc == 3)
            edget_polytope = SEGMENT
        end

        basis, reffe_args, reffe_kwargs = reffe
        pindex_to_cfvertex_to_fvertex = Gridap.ReferenceFEs.get_vertex_permutations(facet_polytope)

        if (Dc == 3)
            edge_reffe = ReferenceFE(edget_polytope, basis, reffe_args...; reffe_kwargs...)
            pindex_to_cevertex_to_evertex = Gridap.ReferenceFEs.get_vertex_permutations(edget_polytope)
        end

        owner_faces_pindex = Vector{Vector{Int}}(undef, Dc - 1)
        owner_faces_lids = Vector{Dict{Int,Tuple{Int,Int,Int}}}(undef, Dc - 1)

        lface_to_cvertices = Gridap.ReferenceFEs.get_faces(Dc == 2 ? QUAD : HEX, Dc - 1, 0)
        owner_faces_pindex[Dc-1], owner_faces_lids[Dc-1] = _compute_owner_faces_pindex_and_lids(Dc,
            num_hanging_faces[Dc],
            hanging_faces_glue[Dc],
            hanging_faces_to_cell[Dc],
            hanging_faces_to_lface[Dc],
            gridap_cell_faces[1],
            gridap_cell_faces[Dc],
            lface_to_cvertices,
            pindex_to_cfvertex_to_fvertex)

        if (Dc == 3)
          owner_faces_pindex[1], owner_faces_lids[1]=
            _compute_owner_edges_pindex_and_lids(
                num_hanging_faces[2],
                hanging_faces_glue[2],
                hanging_faces_to_cell[2],
                hanging_faces_to_lface[2],
                gridap_cell_faces[1],
                gridap_cell_faces[2])
        end
     
        face_dofs_permutations = Vector{Vector{Vector{Int}}}(undef, Dc-1)
        face_dofs_permutations[Dc-1] = 
            first(get_face_dofs_permutations(cell_reffe, Dc-1))
        if (Dc == 3)
           face_dofs_permutations[1] = 
             first(get_face_dofs_permutations(cell_reffe, 1))
        end                

        subface_own_dofs = Vector{Vector{Vector{Int}}}(undef, Dc - 1)
        subface_own_dofs[Dc-1] = _restrict_face_dofs_to_face_dim(cell_reffe,Dc-1)
        if (Dc == 3)
            subface_own_dofs[1] = _restrict_face_dofs_to_face_dim(cell_reffe,1)
        end
        if (_num_face_own_dofs(cell_reffe,0)>0)
            _generate_constraints!(0,
                Dc,
                [gridap_cell_faces[i] for i = 1:Dc],
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
                face_dofs_permutations,
                owner_faces_pindex,
                owner_faces_lids,
                ref_constraints,
                sDOF_to_dof,
                sDOF_to_dofs,
                sDOF_to_coeffs)
        end

        if (Dc == 3)
            if (_num_face_own_dofs(cell_reffe,1)>0)
                _generate_constraints!(1,
                    Dc,
                    [gridap_cell_faces[i] for i = 1:Dc],
                    num_hanging_faces[2],
                    hanging_faces_to_cell[2],
                    hanging_faces_to_lface[2],
                    hanging_faces_owner_face_dofs[2],
                    hanging_faces_glue[2],
                    face_subface_ldof_to_cell_ldof,
                    face_dofs,
                    face_own_dofs,
                    subface_own_dofs,
                    cell_dof_ids,
                    face_dofs_permutations,
                    owner_faces_pindex,
                    owner_faces_lids,
                    ref_constraints,
                    sDOF_to_dof,
                    sDOF_to_dofs,
                    sDOF_to_coeffs)
            end
        end
        if (_num_face_own_dofs(cell_reffe,Dc-1)>0)
            _generate_constraints!(Dc - 1,
                Dc,
                [gridap_cell_faces[i] for i = 1:Dc],
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
                face_dofs_permutations,
                owner_faces_pindex,
                owner_faces_lids,
                ref_constraints,
                sDOF_to_dof,
                sDOF_to_dofs,
                sDOF_to_coeffs)
        end
        sDOF_to_dof, Gridap.Arrays.Table(sDOF_to_dofs), Gridap.Arrays.Table(sDOF_to_coeffs)
    end |> tuple_of_arrays
end

# An auxiliary function which we use in order to generate a version of  
# get_cell_dof_ids() for FE spaces with linear constraints which is suitable 
# for the algorithm which generates the global DoFs identifiers
function fe_space_with_linear_constraints_cell_dof_ids(Uc::FESpaceWithLinearConstraints,sDOF_to_dof)
    U_cell_dof_ids = Gridap.Arrays.Table(get_cell_dof_ids(Uc.space))
    ndata = U_cell_dof_ids.ptrs[end] - 1
    Uc_cell_dof_ids_data = zeros(eltype(U_cell_dof_ids.data), ndata)
    max_negative_minus_one = -maximum(-U_cell_dof_ids.data) - 1
    # max_negative_minus_one can only be zero whenever there are no 
    # negative values in U_cell_dof_ids.data (i.e., no Dirichlet DoFs)
    if (max_negative_minus_one==0) 
        max_negative_minus_one = -1
    end 
    n_cells = length(U_cell_dof_ids)
    n_fdofs = num_free_dofs(Uc.space)
    n_fmdofs = Uc.n_fmdofs
    dof_to_sDOF = Dict(val=>key for (key,val) in enumerate(sDOF_to_dof))
    for cell in 1:n_cells
        pini = U_cell_dof_ids.ptrs[cell]
        pend = U_cell_dof_ids.ptrs[cell+1] - 1
        for p in pini:pend
            dof = U_cell_dof_ids.data[p]
            DOF = Gridap.FESpaces._dof_to_DOF(dof, n_fdofs)
            qini = Uc.DOF_to_mDOFs.ptrs[DOF]
            qend = Uc.DOF_to_mDOFs.ptrs[DOF+1]-1 
            if (!haskey(dof_to_sDOF,dof)) # master DoF
                @assert qend-qini==0
                mDOF = Uc.DOF_to_mDOFs.data[qini]
                mdof = Gridap.FESpaces._DOF_to_dof(mDOF, n_fmdofs)
                Uc_cell_dof_ids_data[p] = mdof
            else # slave DoF
                Uc_cell_dof_ids_data[p] = max_negative_minus_one
            end
        end
    end
    Gridap.Arrays.Table(Uc_cell_dof_ids_data, U_cell_dof_ids.ptrs)
end

function _is_conforming(model::OctreeDistributedDiscreteModel)
    is_local_conforming=map(model.non_conforming_glue) do ncglue 
         all(x->x==0, ncglue.num_hanging_faces)
    end 
    reduction(&,is_local_conforming,init=true,destination=:all).item_ref[]
end 

function _add_constraints(model::GridapDistributed.DistributedDiscreteModel{Dc},
                          reffe,
                          spaces_wo_constraints;
                          conformity=nothing,
                          kwargs...) where {Dc}

    if (_is_conforming(model) || conformity==:L2 )
        spaces_w_constraints=spaces_wo_constraints
        local_cell_dof_ids=map(get_cell_dof_ids,spaces_w_constraints)
    else 
        @assert conformity==nothing || conformity!=:L2
        ref_constraints = _build_constraint_coefficients_matrix_in_ref_space(Dc, reffe)
        face_subface_ldof_to_cell_ldof = Vector{Vector{Vector{Vector{Int32}}}}(undef, Dc-1)
        face_subface_ldof_to_cell_ldof[Dc-1] = _generate_face_subface_ldof_to_cell_ldof(Dc-1, Dc, reffe)
        if (Dc == 3)
            face_subface_ldof_to_cell_ldof[1] =
                _generate_face_subface_ldof_to_cell_ldof(1, Dc, reffe)
        end
        sDOF_to_dof, sDOF_to_dofs, sDOF_to_coeffs =
            generate_constraints(model, spaces_wo_constraints, reffe, ref_constraints, face_subface_ldof_to_cell_ldof)
        spaces_w_constraints = map(spaces_wo_constraints,
            sDOF_to_dof,
            sDOF_to_dofs,
            sDOF_to_coeffs) do V, sDOF_to_dof, sDOF_to_dofs, sDOF_to_coeffs
            @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))]: fe_space_wo_constraints_cell_dof_ids=$(get_cell_dof_ids(V))"
            Vc = FESpaceWithLinearConstraints(sDOF_to_dof, sDOF_to_dofs, sDOF_to_coeffs, V)
        end
        local_cell_dof_ids = map(spaces_w_constraints,sDOF_to_dof) do Vc,sDOF_to_dof
            result = fe_space_with_linear_constraints_cell_dof_ids(Vc,sDOF_to_dof)
            @debug "[$(MPI.Comm_rank(MPI.COMM_WORLD))]: fe_space_with_linear_constraints_cell_dof_ids=$(result)"
            result
        end
    end
    nldofs = map(num_free_dofs,spaces_w_constraints)
    cell_gids = get_cell_gids(model)
    gids = GridapDistributed.generate_gids(cell_gids,local_cell_dof_ids,nldofs)
    map(partition(gids)) do indices 
        @debug "[$(part_id(indices))]: l2g_cell_gids=$(local_to_global(indices))"
        @debug "[$(part_id(indices))]: l2o_owner=$(local_to_owner(indices))"
    end
    trian = Triangulation(model)
    vector_type = GridapDistributed._find_vector_type(spaces_w_constraints,gids)
    GridapDistributed.DistributedSingleFieldFESpace(spaces_w_constraints,gids,trian,vector_type)
end

# Generates a new DistributedSingleFieldFESpace composed 
# by local FE spaces with linear multipoint constraints added
function Gridap.FESpaces.FESpace(model::OctreeDistributedDiscreteModel{Dc}, 
                                 reffe::Tuple{Gridap.ReferenceFEs.Lagrangian,Any,Any}; 
                                 kwargs...) where {Dc}
    spaces_wo_constraints = map(local_views(model)) do m
        FESpace(m, reffe; kwargs...)
    end
    _add_constraints(model,reffe,spaces_wo_constraints;kwargs...)
end

function Gridap.FESpaces.FESpace(model::OctreeDistributedDiscreteModel{Dc}, 
                                 reffe::Tuple{Gridap.ReferenceFEs.RaviartThomas,Any,Any}; 
                                 conformity=nothing,kwargs...) where {Dc}

    cell_reffes = map(local_views(model.dmodel)) do m
        basis,reffe_args,reffe_kwargs = reffe
        cell_reffe = ReferenceFE(m,basis,reffe_args...;reffe_kwargs...)
    end
    sign_flips=GridapDistributed._generate_sign_flips(model.dmodel,cell_reffes)
    spaces_wo_constraints = map(local_views(model.dmodel),sign_flips,cell_reffes) do m,sign_flip,cell_reffe
       conf = Conformity(Gridap.Fields.testitem(cell_reffe),conformity)
       cell_fe = CellFE(m,cell_reffe,conf,sign_flip)
       FESpace(m, cell_fe; kwargs...)
    end
    _add_constraints(model,reffe,spaces_wo_constraints;conformity=conformity,kwargs...)
end
