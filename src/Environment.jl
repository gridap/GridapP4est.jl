# This variable is necessary to avoid libsc being initialized twice
const _INITIALIZED = Ref(false)

"""
See P4est_wrapper.jl/src/bindings/sc_common.jl for possible/valid
argument values for the p4est_verbosity_level parameter
"""
function Init(parts::MPIArray;p4est_verbosity_level=P4est_wrapper.SC_LP_DEFAULT)
  if !MPI.Initialized()
    @error "MPI not Initialized!"
  end
  if _INITIALIZED[] == false
    sc_init(parts.comm, Cint(false), Cint(false), C_NULL, p4est_verbosity_level)
    p4est_init(C_NULL, p4est_verbosity_level)
    _INITIALIZED[] = true
  end
  return nothing
end

function Finalize()
  GC.gc() # Finalize all objects out of scope at this point

  # The function call to sc_finalize() is useful as, among others, it also
  # double checks for memory leaks as a result of improper usage
  # of p4est/sc; if there are leaks, then it generates a fatal
  # error.

  # HOWEVER ... we cannot call it as we have realized there are
  # scenarios where the call to GC.gc() above does not necessarily
  # run synchronously, but scheduled in a co-routine/task

  # As a result, sc_finalize() may generate an error and is not safe to
  # to call it at this point.

  # sc_finalize()
  return nothing
end

function with(f,parts::MPIArray;kwargs...)
  Init(parts;kwargs...)
  out = f()
  Finalize()
  return out
end

function num_parts(comm::MPI.Comm)
  if comm != MPI.COMM_NULL
    nparts = MPI.Comm_size(comm)
  else
    nparts = -1
  end
  nparts
end

function get_part_id(comm::MPI.Comm)
  if comm != MPI.COMM_NULL
    id = MPI.Comm_rank(comm)+1
  else
    id = -1
  end
  id
end

function i_am_in(comm::MPI.Comm)
  get_part_id(comm) >=0
end

function i_am_in(comm::MPIArray)
  i_am_in(comm.comm)
end

# This type is required because MPIArray from PArrays 
# cannot be instantiated with a NULL communicator
struct MPIVoidVector{T} <: AbstractVector{T}
    comm::MPI.Comm
    function MPIVoidVector(::Type{T}) where {T}
        new{T}(MPI.COMM_NULL)
    end
end

Base.size(a::MPIVoidVector) = (0,)
Base.IndexStyle(::Type{<:MPIVoidVector}) = IndexLinear()
function Base.getindex(a::MPIVoidVector,i::Int)
  error("Indexing of MPIVoidVector not possible.")
end
function Base.setindex!(a::MPIVoidVector,v,i::Int)
  error("Indexing of MPIVoidVector not possible.")
end
function Base.show(io::IO,k::MIME"text/plain",data::MPIVoidVector)
  println(io,"MPIVoidVector")
end

function i_am_in(comm::MPIVoidVector)
  i_am_in(comm.comm)
end