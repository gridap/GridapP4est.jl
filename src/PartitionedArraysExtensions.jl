
function PartitionedArrays.with_backend(driver::Function,b::MPIBackend,nparts::Union{Int,NTuple{N,Int} where N},args...;kwargs...)
  if !MPI.Initialized()
    MPI.Init()
  end
  if MPI.Comm_size(MPI.COMM_WORLD) == 1
    part = get_part_ids(b,nparts)
    driver(part,args...;kwargs...)
  else
    try
       part = get_part_ids(b,nparts)
       if i_am_in(part) 
         driver(part,args...;kwargs...)
       end
    catch e
      @error "" exception=(e, catch_backtrace())
      if MPI.Initialized() && !MPI.Finalized()
        MPI.Abort(MPI.COMM_WORLD,1)
      end
    end
  end
  # We are NOT invoking MPI.Finalize() here because we rely on
  # MPI.jl, which registers MPI.Finalize() in atexit()
end

function generate_level_parts(parts,num_procs_x_level)
  root_comm = parts.comm
  rank = MPI.Comm_rank(root_comm)
  size = MPI.Comm_size(root_comm)
  Gridap.Helpers.@check all(num_procs_x_level .<= size)
  Gridap.Helpers.@check all(num_procs_x_level .>= 1)

  @static if isdefined(MPI,:MPI_UNDEFINED)
    mpi_undefined = MPI.MPI_UNDEFINED[]
  else
    mpi_undefined = MPI.API.MPI_UNDEFINED[]
  end
  
  nlevs = length(num_procs_x_level)
  level_parts=Vector{typeof(parts)}(undef,nlevs)
  for l = 1:nlevs
    lsize = num_procs_x_level[l]
    if (l > 1) && (lsize == num_procs_x_level[l-1])
      level_parts[l] = level_parts[l-1]
    else
      if rank < lsize
        comm=MPI.Comm_split(root_comm,0,0)
      else
        comm=MPI.Comm_split(root_comm,mpi_undefined,mpi_undefined)
      end
      level_parts[l] = get_part_ids(comm)
    end
  end
  return level_parts
end
