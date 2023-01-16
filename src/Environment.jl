
const _NREFS = Ref(0)

"""
See P4est_wrapper.jl/src/bindings/sc_common.jl for possible/valid
argument values for the p4est_verbosity_level parameter
"""
function Init(parts::MPIData;p4est_verbosity_level=P4est_wrapper.SC_LP_DEFAULT)
  if !MPI.Initialized()
    @error "MPI not Initialized!"
  end

  sc_init(parts.comm, Cint(true), Cint(true), C_NULL, p4est_verbosity_level)
  p4est_init(C_NULL, p4est_verbosity_level)

  return nothing
end

function Finalize()
  GC.gc() # Finalize all object out of scope at this point
  if _NREFS[] != 0
    @warn "$(_NREFS[]) object(s) still not finalized before calling GridapP4est.Finalize()"
  end
  _NREFS[] = 0
  # This function call is useful, as among others, it also double checks
  # for memory leaks as a result of improper usage of p4est/sc. If there are leaks, 
  # then it generates a fatal error.
  sc_finalize()
  return nothing
end

function with(f,parts::MPIData;kwargs...)
  Init(parts;kwargs...)
  out = f()
  Finalize()
  return out
end