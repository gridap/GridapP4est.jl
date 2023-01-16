"""
See P4est_wrapper.jl/src/bindings/sc_common.jl for possible/valid
argument values for the p4est_verbosity_level parameter
"""
function Init(parts::MPIData;p4est_verbosity_level=P4est_wrapper.SC_LP_DEFAULT)
  if !MPI.Initialized()
    @error "MPI not Initialized!"
  end
  sc_init(parts.comm, Cint(false), Cint(false), C_NULL, p4est_verbosity_level)
  p4est_init(C_NULL, p4est_verbosity_level)
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

function with(f,parts::MPIData;kwargs...)
  Init(parts;kwargs...)
  out = f()
  Finalize()
  return out
end
