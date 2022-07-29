


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


function GMG!(x::PVector,
              b::PVector,
              mh::ModelHierarchy,
              smatrices,
              interpolations,
              restrictions;
              rtol=1.0e-06,
              maxiter=100,
              pre_smoothers::AbstractVector{<:Gridap.Algebra.LinearSolver}=Fill(RichardsonSmoother(JacobiLinearSolver(),10),num_levels(mh)-1),
              post_smoothers::AbstractVector{<:Gridap.Algebra.LinearSolver}=pre_smoothers,
              coarsest_solver::Gridap.Algebra.LinearSolver=BackslashSolver(),
              verbose=false)

  work_vectors=allocate_work_vectors(mh,smatrices)
  pre_smoothers_caches=setup_smoothers_caches(mh,pre_smoothers,smatrices)
  if (!(pre_smoothers===post_smoothers))
    post_smoothers_caches=setup_smoothers_caches(mh,post_smoothers,smatrices)
  else
    post_smoothers_caches=pre_smoothers_caches
  end
  coarsest_solver_cache=setup_coarsest_solver_cache(mh,coarsest_solver,smatrices)

  Ah=smatrices[1]
  rh = PVector(0.0,Ah.rows)
  rh .= b .- Ah*x
  nrm_r0 = norm(rh)
  nrm_r  = nrm_r0
  current_iter = 0
  rel_res = nrm_r / nrm_r0
  model = get_level_model(mh,1)

  if (GridapP4est.i_am_main(model.parts))
    @printf "%6s  %12s" "Iter" "Rel res\n"
    @printf "%6i  %12.4e\n" current_iter rel_res
  end

  while current_iter <= maxiter && rel_res > rtol
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
