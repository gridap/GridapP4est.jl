using MPI
using PartitionedArrays
using Gridap
using GridapP4est
using P4est_wrapper

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



function run(parts,num_procs_x_level)
  root_comm = parts.comm
  rank = MPI.Comm_rank(root_comm)
  level_parts = generate_level_parts(parts,num_procs_x_level)

  domain=(0,1,0,1)
  subdomains=(1,1)
  num_uniform_refinements=1
  coarse_discrete_model=CartesianDiscreteModel(domain,subdomains)
  model=OctreeDistributedDiscreteModel(level_parts[1],
                                       coarse_discrete_model)
  model,glue=refine(model)
  model,glue=refine(model)
  model=redistribute(model,level_parts[2])
end

# Give me how many processors you want per level
# in an array with as many entries as levels
# num_procs_x_level = [4,4,4,2,2,2,1]
num_procs_x_level = [1,2]
if !MPI.Initialized()
  MPI.Init()
end
parts = get_part_ids(mpi,2)
run(parts,num_procs_x_level)
MPI.Finalize()
