
struct GMGLinearSolver{A,B,C,D,E,F,G,H} <: Gridap.Algebra.LinearSolver
  mh              :: ModelHierarchy
  smatrices       :: A
  interp          :: B
  restrict        :: C
  pre_smoothers   :: D
  post_smoothers  :: E
  coarsest_solver :: F
  maxiter         :: G
  rtol            :: H
  verbose         :: Bool
  mode            :: Symbol
end

function GMGLinearSolver(mh,
      smatrices,
      interp,
      restrict;
      pre_smoothers=Fill(RichardsonSmoother(JacobiLinearSolver(),10),num_levels(mh)-1),
      post_smoothers=pre_smoothers,
      coarsest_solver=Gridap.Algebra.BackslashSolver(),
      maxiter=100,
      rtol=1.0e-06,
      verbose::Bool=false,
      mode=:preconditioner)

  Gridap.Helpers.@check mode == :preconditioner || mode == :solver
  Gridap.Helpers.@check isa(maxiter,Integer)
  Gridap.Helpers.@check isa(rtol,Real)

  A=typeof(smatrices)
  B=typeof(interp)
  C=typeof(restrict)
  D=typeof(pre_smoothers)
  E=typeof(post_smoothers)
  F=typeof(coarsest_solver)
  G=typeof(maxiter)
  H=typeof(rtol)
  GMGLinearSolver{A,B,C,D,E,F,G,H}(mh,
                   smatrices,
                   interp,
                   restrict,
                   pre_smoothers,
                   post_smoothers,
                   coarsest_solver,
                   maxiter,
                   rtol,
                   verbose,
                   mode)
end

struct GMGSymbolicSetup <: Gridap.Algebra.SymbolicSetup
  solver :: GMGLinearSolver
end

function Gridap.Algebra.symbolic_setup(solver::GMGLinearSolver,mat::AbstractMatrix)
  GMGSymbolicSetup(solver)
end

struct GMGNumericalSetup{A,B,C,D} <: Gridap.Algebra.NumericalSetup
  solver                 :: GMGLinearSolver
  pre_smoothers_caches   :: A
  post_smoothers_caches  :: B
  coarsest_solver_cache  :: C
  work_vectors           :: D

  function GMGNumericalSetup(ss::GMGSymbolicSetup)
    mh              = ss.solver.mh
    pre_smoothers   = ss.solver.pre_smoothers
    post_smoothers  = ss.solver.post_smoothers
    smatrices       = ss.solver.smatrices
    coarsest_solver = ss.solver.coarsest_solver

    work_vectors=allocate_work_vectors(mh,smatrices)
    pre_smoothers_caches=setup_smoothers_caches(mh,pre_smoothers,smatrices)
    if (!(pre_smoothers===post_smoothers))
      post_smoothers_caches=setup_smoothers_caches(mh,post_smoothers,smatrices)
    else
      post_smoothers_caches=pre_smoothers_caches
    end
    coarsest_solver_cache=setup_coarsest_solver_cache(mh,coarsest_solver,smatrices)
    A=typeof(pre_smoothers_caches)
    B=typeof(post_smoothers_caches)
    C=typeof(coarsest_solver_cache)
    D=typeof(work_vectors)
    new{A,B,C,D}(ss.solver,
                 pre_smoothers_caches,
                 post_smoothers_caches,
                 coarsest_solver_cache,
                 work_vectors)
  end
end

function Gridap.Algebra.numerical_setup(ss::GMGSymbolicSetup,mat::AbstractMatrix)
  GMGNumericalSetup(ss)
end

function setup_smoothers_caches(mh,smoothers,smatrices)
  Gridap.Helpers.@check length(smoothers) == num_levels(mh)-1
  nlevs=num_levels(mh)
  # Last (i.e., coarsest) level does not need pre-/post-smoothing
  caches=Vector{Any}(undef,nlevs-1)
  for i=1:nlevs-1
    model = get_level_model(mh,i)
    if (GridapP4est.i_am_in(model.parts))
      ss = symbolic_setup(smoothers[i], smatrices[i])
      caches[i] = numerical_setup(ss, smatrices[i])
    end
  end
  caches
end

function setup_coarsest_solver_cache(mh,coarsest_solver,smatrices)
  cache=nothing
  nlevs=num_levels(mh)
  model=get_level_model(mh,nlevs)
  if (GridapP4est.i_am_in(model.parts))
    if (num_parts(model.parts)==1)
      cache=map_parts(smatrices[nlevs].owned_owned_values) do Ah
        ss = symbolic_setup(coarsest_solver, Ah)
        numerical_setup(ss, Ah)
      end
      cache=cache.part
    else
      ss = symbolic_setup(coarsest_solver, smatrices[nlevs])
      cache=numerical_setup(ss, smatrices[nlevs])
    end
  end
  cache
end

function allocate_level_work_vectors(mh,smatrices,lev)
  modelH = get_level_model(mh,lev+1)
  dxh  = PVector(0.0, smatrices[lev].cols)
  Adxh = PVector(0.0, smatrices[lev].rows)
  rh   = PVector(0.0, smatrices[lev].rows)
  if (GridapP4est.i_am_in(modelH.parts))
    AH  = smatrices[lev+1]
    rH  = PVector(0.0,AH.cols)
    dxH = PVector(0.0,AH.cols)
  else
    rH  = nothing
    dxH = nothing
  end
  (dxh,Adxh,dxH,rH)
end

function allocate_work_vectors(mh,smatrices)
  nlevs=num_levels(mh)
  work_vectors=Vector{Any}(undef,nlevs-1)
  for i=1:nlevs-1
    model = get_level_model(mh,i)
    if (GridapP4est.i_am_in(model.parts))
      work_vectors[i]=allocate_level_work_vectors(mh,smatrices,i)
    end
  end
  work_vectors
end

function apply_GMG_level!(xh,
                          rh,
                          lev,
                          mh,
                          smatrices,
                          restrictions,
                          interpolations,
                          pre_smoothers_caches,
                          post_smoothers_caches,
                          coarsest_solver_cache,
                          work_vectors;
                          verbose=false)

  modelh = get_level_model(mh,lev)
  if (GridapP4est.i_am_in(modelh.parts))
    if (lev==num_levels(mh))
      if (GridapP4est.num_parts(modelh.parts)==1)
         map_parts(smatrices[lev].owned_owned_values,
                   xh.owned_values,
                   rh.owned_values) do Ah, xh, rh
            solve!(xh,coarsest_solver_cache,rh)
            rh .= rh - Ah*xh
         end
      else
        solve!(xh,coarsest_solver_cache,rh)
      end
    else
      Ah = smatrices[lev]
      (dxh,Adxh,dxH,rH)=work_vectors[lev]

      # Pre-smooth current solution
      solve!(xh, pre_smoothers_caches[lev], rh)

      # Restrict the residual
      mul!(rH,restrictions[lev],rh; verbose=verbose)

      if dxH!=nothing
        fill!(dxH,0.0)
      end

      # Apply next_level
      apply_GMG_level!(dxH,
                      rH,
                      lev+1,
                      mh,
                      smatrices,
                      restrictions,
                      interpolations,
                      pre_smoothers_caches,
                      post_smoothers_caches,
                      coarsest_solver_cache,
                      work_vectors;
                      verbose=verbose)

      # Interpolate dxH in finer space
      mul!(dxh, interpolations[lev], dxH; verbose=verbose)

      # Update solution
      xh .= xh .+ dxh
      # Update residual
      mul!(Adxh, Ah, dxh)
      rh .= rh .- Adxh

      # Post-smooth current solution
      solve!(xh, post_smoothers_caches[lev], rh)

    end
  end
end

function Gridap.Algebra.solve!(
  x::Vector,ns::GMGNumericalSetup,b::Vector)
  smatrices = ns.solver.smatrices
  Ah        = smatrices[1]
  px        = PVector(0.0, Ah.cols)
  pb        = PVector(0.0, Ah.rows)
  Gridap.Helpers.@check num_parts(pb.values) == 1
  Gridap.Helpers.@check num_parts(px.values) == 1
  px.values.part .= x
  pb.values.part .= b
  Gridap.Algebra.solve!(px,ns,pb)
  x .= px.values.part
end


function Gridap.Algebra.solve!(
  x::AbstractVector,ns::GMGNumericalSetup,b::AbstractVector)

  smatrices      = ns.solver.smatrices
  mh             = ns.solver.mh
  maxiter        = ns.solver.maxiter
  rtol           = ns.solver.rtol
  restrictions   = ns.solver.restrict
  interpolations = ns.solver.interp
  verbose        = ns.solver.verbose
  mode           = ns.solver.mode

  pre_smoothers_caches  = ns.pre_smoothers_caches
  post_smoothers_caches = ns.post_smoothers_caches
  coarsest_solver_cache = ns.coarsest_solver_cache
  work_vectors          = ns.work_vectors

  Ah=smatrices[1]

  if (mode==:preconditioner)
    fill!(x,0.0)
    rh = copy(b)
  else
    rh = PVector(0.0,Ah.rows)
    rh .= b .- Ah*x
  end

  nrm_r0 = norm(rh)
  nrm_r  = nrm_r0
  current_iter = 0
  rel_res = nrm_r / nrm_r0
  model = get_level_model(mh,1)

  if (GridapP4est.i_am_main(model.parts))
    @printf "%6s  %12s" "Iter" "Rel res\n"
    @printf "%6i  %12.4e\n" current_iter rel_res
  end

  while current_iter < maxiter && rel_res > rtol
    apply_GMG_level!(x,
                     rh,
                     1,
                     mh,
                     smatrices,
                     restrictions,
                     interpolations,
                     pre_smoothers_caches,
                     post_smoothers_caches,
                     coarsest_solver_cache,
                     work_vectors;
                     verbose=verbose)
    nrm_r = norm(rh)
    rel_res = nrm_r / nrm_r0
    current_iter += 1
    if (GridapP4est.i_am_main(model.parts))
      @printf "%6i  %12.4e\n" current_iter rel_res
    end
  end
  converged=(rel_res < rtol)
  return current_iter, converged
end

function LinearAlgebra.ldiv!(x::AbstractVector,ns::GMGNumericalSetup,b::AbstractVector)
  solve!(x,ns,b)
end
