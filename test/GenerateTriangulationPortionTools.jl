function _generate_triangulation_portion(ranks, fmodel)
        trians = map(ranks, 
             local_views(fmodel.dmodel), 
                         partition(get_cell_gids(fmodel))) do rank, lmodel, indices
            mask = Vector{Bool}(undef,num_cells(lmodel))
            mask .= false
            cell_to_part = local_to_owner(indices)
            graph = GridapDistributed.compute_cell_graph(lmodel)
            ncells = num_cells(lmodel)
            icell_to_jcells_ptrs = graph.colptr
            icell_to_jcells_data = graph.rowval
            for icell in 1:ncells
               if cell_to_part[icell] == rank
                  pini = icell_to_jcells_ptrs[icell]
                  pend = icell_to_jcells_ptrs[icell+1]-1
                  for p in pini:pend
                     jcell = icell_to_jcells_data[p]
                     if cell_to_part[jcell] != rank
                        mask[icell] = true
                     end                 
                  end
               end
            end
            Triangulation(lmodel, mask)
        end
        GridapDistributed.DistributedTriangulation(trians,fmodel)
  end

  function _generate_triangulation_portion(ranks, fmodel, ctrian, glue)
        trians = map(ranks,
                 local_views(fmodel.dmodel),
                 local_views(ctrian),
                 glue) do rank, fmodel, ctrian, glue
            mask = Vector{Bool}(undef,num_cells(fmodel))
            mask .= false

            cglue = get_glue(ctrian, Val{num_cell_dims(ctrian)}())
            for ccell in cglue.tface_to_mface
                fcells = glue.o2n_faces_map[ccell]
                for fcell in fcells
                    mask[fcell] = true
                end
            end
            Triangulation(fmodel, mask)
        end
        GridapDistributed.DistributedTriangulation(trians,fmodel)
  end

  function generate_triangulation_portion(ranks,fmodel; ctrian=nothing, glue=nothing)
    if (length(ranks)==1 && ctrian==nothing)
       trians = map(ranks, local_views(fmodel.dmodel)) do rank, fmodel
            mask = Vector{Bool}(undef,num_cells(fmodel))
            mask .= false
            for cell=1:Int(round(num_cells(fmodel)*0.25))
               mask[cell] = true
            end
            for cell=Int(round(num_cells(fmodel)*0.75)):num_cells(fmodel)
               mask[cell] = true
            end
            Triangulation(fmodel, mask)
        end
        return GridapDistributed.DistributedTriangulation(trians,fmodel)
    else 
        if ctrian==nothing
            @assert glue==nothing
            return _generate_triangulation_portion(ranks, fmodel)
        else
            return _generate_triangulation_portion(ranks, fmodel, ctrian, glue)
        end
    end
  end
