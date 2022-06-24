function PartitionedArrays.num_parts(comm::MPI.Comm)
  if comm != MPI.COMM_NULL
    nparts = MPI.Comm_size(comm)
  else
    nparts = -1
  end
  nparts
end

function PartitionedArrays.get_part_id(comm::MPI.Comm)
  if comm != MPI.COMM_NULL
    id = MPI.Comm_rank(comm)+1
  else
    id = -1
  end
  id
end

function i_am_in(comm::MPI.Comm)
  PartitionedArrays.get_part_id(comm) >=0
end

function PartitionedArrays.get_part_ids(comm::MPI.Comm)
  rank = PartitionedArrays.get_part_id(comm)
  nparts = PartitionedArrays.num_parts(comm)
  PartitionedArrays.MPIData(rank,comm,(nparts,))
end

function generate_level_parts(parts,num_procs_x_level)
  root_comm = parts.comm
  rank = MPI.Comm_rank(root_comm)
  size = MPI.Comm_size(root_comm)
  Gridap.Helpers.@check all(num_procs_x_level .<= size)
  Gridap.Helpers.@check all(num_procs_x_level .>= 1)

  level_parts=Vector{typeof(parts)}(undef,length(num_procs_x_level))
  for l=1:length(num_procs_x_level)
    lsize=num_procs_x_level[l]
    if l>1 && lsize==num_procs_x_level[l-1]
      level_parts[l]=level_parts[l-1]
    else
      if rank < lsize
        comm=MPI.Comm_split(root_comm, 0, 0)
      else
        comm=MPI.Comm_split(root_comm, MPI.MPI_UNDEFINED, MPI.MPI_UNDEFINED)
      end
      level_parts[l]=get_part_ids(comm)
    end
  end
  return level_parts
end
