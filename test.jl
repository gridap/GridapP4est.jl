using Gridap 
using GridapP4est 
using MPI 
using PartitionedArrays
using GridapDistributed
using Logging

debug_logger = ConsoleLogger(stderr, Logging.Debug)
global_logger(debug_logger); # Enable the debug logger globally

MPI.Init()

ranks  = with_mpi() do distribute
  distribute(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))
end

coarse_model=CartesianDiscreteModel((0,1,0,1),(1,1))

model=AnisotropicallyAdapted3DDistributedDiscreteModel(ranks,coarse_model,1,1)

writevtk(model, "model");

ref_coarse_flags=map(ranks,partition(get_cell_gids(model.dmodel))) do rank,indices
  flags=zeros(Int,length(indices))
  flags.=nothing_flag        
  flags[1]=refine_flag
  flags[end]=refine_flag
  flags
end

model,glue=GridapP4est.vertically_adapt(model,ref_coarse_flags);

writevtk(Triangulation(model), "trian");