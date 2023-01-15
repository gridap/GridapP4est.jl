
const _NREFS = Ref(0)

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
    @warn "$(_NREFS[]) objects still not finalized before calling GridapP4est.Finalize()"
  end
  _NREFS[] = 0
  return nothing
end

function with(f,parts::MPIData;kwargs...)
  Init(parts;kwargs...)
  out = f()
  Finalize()
  return out
end