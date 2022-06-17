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
