
const _NREFS = Ref(0)
const _INITIALIZED = Ref(false)

function Init(parts::MPIData;p4est_verbosity_level=P4est_wrapper.SC_LP_DEFAULT,finalize_atexit=true)
  if !MPI.Initialized()
      MPI.Init()
  end

  if finalize_atexit
    atexit(Finalize)
  end
  Finalize()

  # sc_init(parts.comm, Cint(true), Cint(true), C_NULL, p4est_verbosity_level)
  p4est_init(C_NULL, p4est_verbosity_level)
  _INITIALIZED[] = true

  return nothing
end

function Initialized()
  return _INITIALIZED[] == true
end

function Finalize()
  if Initialized()
    GC.gc() # Finalize all object out of scope at this point
    if _NREFS[] != 0
      @warn "$(_NREFS[]) objects still not finalized before calling GridapP4est.Finalize()"
    end
    _NREFS[] = 0
    _INITIALIZED[] = false
    #sc_finalize()
  end
  return nothing
end

function with(f,parts::MPIData;kwargs...)
  Init(parts;kwargs...)
  out = f()
  Finalize()
  return out
end