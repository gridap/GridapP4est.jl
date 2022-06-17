# VERY IMPORTANT NOTE: This was WIP, just to explore DMPlex machinery.
# We have abandoned temporarily this approach. At present, the code written here
# just does not work. I am committing it so that we keep track of it.

module DMPlexTests
  using MPI
  using Gridap
  using Gridap.Algebra
  using Gridap.FESpaces
  using PartitionedArrays
  using GridapDistributed
  using GridapP4est
  using GridapPETSc
  using Test
  using ArgParse
  using P4est_wrapper

  struct DMPlexToGridapDistributed{T}
    petsc_dmplex_face_gids::Vector{T}
  end

  function setup_dmplex(parts::MPIData{<:Integer},model)
    comm=parts.comm
    dm=Ref{GridapP4est.DM}()
    @check_error_code GridapP4est.DMPlexCreate(comm,dm)

    Dc=num_cell_dims(model)
    nofaces=0
    for d=1:Dc+1
      face_gids=get_face_gids(model,d-1)
      nodfaces = map_parts(num_oids,face_gids.partition)
      map_parts(nodfaces) do nodfaces
        nofaces=nofaces+nodfaces
        nothing
      end
    end
    nofaces=map_parts(parts) do part
      println("$(part) $(nofaces)")
      nofaces
    end
    first_gface,ngfaces= xscan(+,reduce,nofaces,init=0)

    current_gface=0
    map_parts(first_gface) do first_gface
      current_gface=first_gface
    end

    petsc_dmplex_face_gids=Vector{PVector{Int}}(undef,Dc+1)
    for d=1:Dc+1
      face_gids=get_face_gids(model,d-1)
      petsc_dmplex_face_gids[d]=PVector{Int}(undef,face_gids)
      map_parts(parts,petsc_dmplex_face_gids[d].values,face_gids.partition) do part, values, gids
        # if (part==2)
        #   println("XXX $(part) $(current_gface)")
        # end
        for i=1:length(gids.lid_to_part)
          if gids.lid_to_part[i]==part
            values[i]=current_gface
            current_gface=current_gface+1
          end
        end
      end
      exchange!(petsc_dmplex_face_gids[d])
    end

    map_parts(first_gface,nofaces) do first_gface, nofaces
       #@check_error_code GridapP4est.DMPlexSetChart(dm[], first_gface, first_gface+nofaces)
       @check_error_code GridapP4est.DMPlexSetChart(dm[], 0, ngfaces)
    end

    for d=2:Dc+1
      faces_gids=get_face_gids(model,d-1)
      map_parts(parts, model.models, faces_gids.partition, petsc_dmplex_face_gids[d].values) do part, model, gids, petsc_dmplex_gids
        topo=Gridap.Geometry.get_grid_topology(model)
        faces=Gridap.Geometry.get_faces(topo,d-1,d-2)
        # if (part==1)
        #   println(petsc_dmplex_gids)
        # end
        for i=1:length(faces.ptrs)-1
          n=faces.ptrs[i+1]-faces.ptrs[i]
          if gids.lid_to_part[i]==part
            fgid=petsc_dmplex_gids[i]
            if (part==1)
               println("$(d) $(PetscInt(fgid)) $(PetscInt(n))")
            end
            @check_error_code GridapP4est.DMPlexSetConeSize(dm[], PetscInt(fgid), PetscInt(n))
          end
        end
      end
    end
    @check_error_code GridapP4est.DMSetUp(dm[])

    carray=Gridap.Arrays.CachedArray(PetscInt[])
    for d=2:Dc+1
     faces_gids=get_face_gids(model,d-1)
     map_parts(parts,
       model.models,
       faces_gids.partition,
       petsc_dmplex_face_gids[d].values,
       petsc_dmplex_face_gids[d-1].values) do part, model, gids, petsc_dmplex_gids_d, petsc_dmplex_gids_dm1
        topo=Gridap.Geometry.get_grid_topology(model)
        faces=Gridap.Geometry.get_faces(topo,d-1,d-2)
        for i=1:length(faces.ptrs)-1
          n=faces.ptrs[i+1]-faces.ptrs[i]
          current=1
          if gids.lid_to_part[i]==part
            fgid=petsc_dmplex_gids_d[i]
            Gridap.Arrays.setsize!(carray,(n,))
            a=carray.array
            for j=faces.ptrs[i]:faces.ptrs[i+1]-1
              flid=faces.data[j]
              a[current]=petsc_dmplex_gids_dm1[flid]
              current=current+1
            end
            if (part==1)
              println("YYY $(d) $(PetscInt(fgid)) $(PetscInt(n)) $(a)")
            end
            @check_error_code GridapP4est.DMPlexSetCone(dm[], PetscInt(fgid), a)
          end
        end
      end
    end
    @check_error_code GridapP4est.DMPlexSymmetrize(dm[])
    @check_error_code GridapP4est.DMPlexStratify(dm[])
    #@check_error_code GridapP4est.DMPlexCheck(dm[])
    #@check_error_code GridapP4est.DMView(dm[],C_NULL)
    #@check_error_code GridapP4est.DMDestroy(dm,C_NULL)
    dm, DMPlexToGridapDistributed(petsc_dmplex_face_gids)
  end

  function setup_petsc_section(parts::MPIData{<:Integer},model,dm,dmplextogdist)
    comm=parts.comm
    sec=Ref{GridapP4est.PetscSection}()

    pStart=Ref{PetscInt}()
    pEnd=Ref{PetscInt}()

    cStart=Ref{PetscInt}()
    cEnd=Ref{PetscInt}()

    eStart=Ref{PetscInt}()
    eEnd=Ref{PetscInt}()

    vStart=Ref{PetscInt}()
    vEnd=Ref{PetscInt}()

    @check_error_code GridapP4est.PetscSectionCreate(comm,sec)
    @check_error_code GridapP4est.DMPlexGetChart(dm[],pStart,pEnd)

    @check_error_code GridapP4est.DMPlexGetHeightStratum(dm[],0,cStart,cEnd)
    @check_error_code GridapP4est.DMPlexGetHeightStratum(dm[],1,eStart,eEnd)
    @check_error_code GridapP4est.DMPlexGetHeightStratum(dm[],2,vStart,vEnd)

    map_parts(parts) do part
      if (part==1)
        println("$(pStart[]) $(pEnd[])")
        println("$(cStart[]) $(cEnd[])")
        println("$(eStart[]) $(eEnd[])")
        println("$(vStart[]) $(vEnd[])")
      end
    end


    @check_error_code GridapP4est.PetscSectionSetChart(sec[], pStart[], pEnd[])
    for c=cStart[]:cEnd[]-1
      @check_error_code GridapP4est.PetscSectionSetDof(sec[], c, 0)
    end
    for v=vStart[]:vEnd[]-1
      @check_error_code GridapP4est.PetscSectionSetDof(sec[], v, 1);
    end
    for e=eStart[]:eEnd[]-1
      @check_error_code GridapP4est.PetscSectionSetDof(sec[], e, 0);
    end
    @check_error_code GridapP4est.PetscSectionSetUp(sec[])
    sec
  end

  function run(parts,subdomains,num_uniform_refinements)
    if length(subdomains)==2
      domain=(0,1,0,1)
    else
      @assert length(subdomains)==3
      domain=(0,1,0,1,0,1)
    end

    coarse_discrete_model=CartesianDiscreteModel(domain,(1,1))
    model=UniformlyRefinedForestOfOctreesDiscreteModel(parts,
                                                       coarse_discrete_model,
                                                       num_uniform_refinements)

    GridapPETSc.with() do
       dm,dmplex2gdist=setup_dmplex(parts,model)
       setup_petsc_section(parts,model,dm,dmplex2gdist)
    end
  end

  if !MPI.Initialized()
    MPI.Init()
  end
  subdomains=(2,2)
  parts = get_part_ids(mpi,(prod(subdomains)))
  run(parts,subdomains,1)
  MPI.Finalize()
end # module
