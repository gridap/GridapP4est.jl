
function _allocate_cell_wise_dofs(cell_to_ldofs)
  map_parts(cell_to_ldofs) do cell_to_ldofs
    cache = array_cache(cell_to_ldofs)
    ncells = length(cell_to_ldofs)
    ptrs = Vector{Int32}(undef,ncells+1)
    for cell in 1:ncells
      ldofs = getindex!(cache,cell_to_ldofs,cell)
      ptrs[cell+1] = length(ldofs)
    end
    PArrays.length_to_ptrs!(ptrs)
    ndata = ptrs[end]-1
    data = Vector{Float64}(undef,ndata)
    PArrays.Table(data,ptrs)
  end
end

function _update_cell_dof_values_with_local_info!(cell_dof_values_new,
                                                  cell_dof_values_old,
                                                  new2old)
   map_parts(cell_dof_values_new,
             cell_dof_values_old,
             new2old) do cell_dof_values_new,cell_dof_values_old,new2old
    ocache=array_cache(cell_dof_values_old)
    for (ncell,ocell) in enumerate(new2old)
      if ocell!=0
        # Copy ocell to ncell
        oentry=getindex!(ocache,cell_dof_values_old,ocell)
        range=cell_dof_values_new.ptrs[ncell]:cell_dof_values_new.ptrs[ncell+1]-1
        cell_dof_values_new.data[range] .= oentry
        # println("YYY $(ncell) $(ocell) $(range) $(cell_dof_values_new.data[range])")
      end
    end
   end
end

function allocate_comm_data(cell_dof_values,lids)
  map_parts(cell_dof_values,lids) do cell_dof_values,lids
    cache = array_cache(cell_dof_values)
    # println("YYY $(length(lids.ptrs))")
    n = length(lids)
    ptrs = Vector{Int32}(undef,n+1)
    ptrs.= 0
    for i=1:n
      for j=lids.ptrs[i]:lids.ptrs[i+1]-1
        # println("XXXX $(n) $(i) $(j) $(length(lids.data))")
        cell=lids.data[j]
        #ldofs=getindex!(cache,cell_dof_values,cell)
        ptrs[i+1]=ptrs[i+1]+4#length(ldofs)
      end
    end
    PArrays.length_to_ptrs!(ptrs)
    ndata = ptrs[end]-1
    data = Vector{Float64}(undef,ndata)
    PArrays.Table(data,ptrs)
  end
end

function pack_snd_data!(snd_data,cell_dof_values,snd_lids)
  map_parts(snd_data,cell_dof_values,snd_lids) do snd_data,cell_dof_values,snd_lids
    cache = array_cache(cell_dof_values)
    s=1
    for i=1:length(snd_lids)
      for j=snd_lids.ptrs[i]:snd_lids.ptrs[i+1]-1
        cell=snd_lids.data[j]
        ldofs=getindex!(cache,cell_dof_values,cell)
        e=s+length(ldofs)-1
        range=s:e
        snd_data.data[range] .= ldofs
        # println("AAA $(cell) $(ldofs) $(range)")
        s=e+1
      end
    end
  end
end

function unpack_rcv_data!(cell_dof_values,rcv_data,rcv_lids)
  map_parts(cell_dof_values,rcv_data,rcv_lids) do cell_dof_values,rcv_data,rcv_lids
    s=1
    #println("AAAA $(rcv_data)")
    for i=1:length(rcv_lids)
      for j=rcv_lids.ptrs[i]:rcv_lids.ptrs[i+1]-1
        cell=rcv_lids.data[j]
        range_cell_dof_values = cell_dof_values.ptrs[cell]:cell_dof_values.ptrs[cell+1]-1
        e=s+length(range_cell_dof_values)-1
        range_rcv_data = s:e
        cell_dof_values.data[range_cell_dof_values] .= rcv_data.data[range_rcv_data]
        #println("BBB $(cell)
        #             $(cell_dof_values.data[range_cell_dof_values])
        #             $(rcv_data.data[range_rcv_data])")
        s=e+1
      end
    end
  end
end

function redistribute_fe_function(uh_old::GridapDistributed.DistributedCellField,
                                  Uh_old::GridapDistributed.DistributedSingleFieldFESpace,
                                  Uh_new::GridapDistributed.DistributedSingleFieldFESpace,
                                  model_new,
                                  glue::RedistributeGlue;
                                  reverse=false)

  if !reverse
    lids_rcv=glue.lids_rcv
    lids_snd=glue.lids_snd
    parts_rcv=glue.parts_rcv
    parts_snd=glue.parts_snd
    new2old=glue.new2old
  else
    lids_rcv=glue.lids_snd
    lids_snd=glue.lids_rcv
    parts_rcv=glue.parts_snd
    parts_snd=glue.parts_rcv
    new2old=glue.old2new
  end

  cell_dof_values_old = map_parts(get_cell_dof_values,uh_old.fields)

  snd_data=allocate_comm_data(cell_dof_values_old,lids_snd)
  rcv_data=allocate_comm_data(cell_dof_values_old,lids_rcv)
  pack_snd_data!(snd_data,cell_dof_values_old,lids_snd)

  # map_parts(snd_data,parts_rcv,parts_snd) do snd_data,parts_rcv,parts_snd
  #  println("YYY $(snd_data) $(parts_rcv) $(parts_snd)")
  # end

  tout= async_exchange!(rcv_data,
                        snd_data,
                        parts_rcv,
                        parts_snd,
                        PArrays._empty_tasks(parts_rcv))

  map_parts(schedule,tout)

  cell_to_ldofs_new   = map_parts(get_cell_dof_ids,Uh_new.spaces)
  cell_dof_values_new = _allocate_cell_wise_dofs(cell_to_ldofs_new)

  # We have to build the owned part of "cell_dof_values_new" out of
  #  1. cell_dof_values_old (for those cells s.t. new2old[:]!=0)
  #  2. cell_dof_values_new_rcv (for those cells s.t. new2old[:]=0)
 _update_cell_dof_values_with_local_info!(cell_dof_values_new,
                                          cell_dof_values_old,
                                          glue.new2old)

 map_parts(wait,tout)

 unpack_rcv_data!(cell_dof_values_new,rcv_data,glue.lids_rcv)

 # map_parts(cell_dof_values_new) do cdvn
 # println("333333333333333 $(cdvn)")
 # end

 fgids=get_cell_gids(model_new)
 exchange!(cell_dof_values_new,fgids.exchanger)

 # map_parts(cell_dof_values_new) do cdvn
 #  println("444 $(cdvn)")
 # end

 free_values = map_parts(cell_dof_values_new,Uh_new.spaces) do cdvn, fspace
   fv=Gridap.FESpaces.gather_free_values(fspace,cdvn)
   # println("111 $(fv)")
   fv
 end
 gids = Uh_new.gids
 free_values=PVector(free_values,gids)
 dirichlet_values = map_parts(cell_dof_values_new,Uh_new.spaces) do cdvn, fspace
  dv=Gridap.FESpaces.gather_dirichlet_values(fspace,cdvn)
  # println("222 $(dv)")
  dv
 end
 FEFunction(Uh_new,free_values,dirichlet_values)
end
