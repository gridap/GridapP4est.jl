
function smooth_cache(A::PSparseMatrix)
  Adx = PVector(0.0,A.rows)
  dx  = PVector(0.0,A.cols)
  inv_diag=map_parts(A.owned_owned_values) do a
    1.0 ./ diag(a)
  end
  (Adx,dx,inv_diag)
end

function smooth!(x::PVector, r::PVector, A::PSparseMatrix; maxiter=1)
  cache=smooth_cache(A)
  smooth!(cache, x, r, A; maxiter=maxiter)
end

# By now, just simple Jacobi smoother !
# We may implement in the future block Gauss Seidel, etc.
# We assume that the residual is up-to-date with x.
function smooth!(cache, x::PVector, r::PVector, A::PSparseMatrix; maxiter=1)
  Adx,dx,inv_diag=cache
  iter=0
  while iter <= maxiter
    map_parts(inv_diag,dx.owned_values,r.owned_values) do inv_diag, dx, r
       dx .= inv_diag .* r
    end
    x .= x .+ dx
    mul!(Adx, A, dx)
    r .= r .- Adx
    iter += 1
  end
end

function setup_smooth_caches(mh,smatrices)
  nlevs=num_levels(mh)
  # Last (i.e., coarsest) level does not need
  # pre-/post-smoothing
  caches=Vector{Any}(undef,nlevs-1)
  for i=1:nlevs-1
    model = get_level_model(mh,i)
    if (GridapP4est.i_am_in(model.parts))
      caches[i]=smooth_cache(smatrices[i])
    end
  end
  caches
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
                          smooth_caches,
                          work_vectors;
                          smooth_iter=2)

  modelh = get_level_model(mh,lev)
  if (GridapP4est.i_am_in(modelh.parts))
    if (lev==num_levels(mh))
      # TO-DO: replace this coarsest solver by a call to a
      #        a sparse direct solver, e.g., PETSc+MUMPs
      # TO-DO: avoid re-computing numerical setup each time
      map_parts(smatrices[lev].owned_owned_values,
                xh.owned_values,
                rh.owned_values) do Ah, xh, rh
         xh .= Ah \ rh
      end
    else
      Ah = smatrices[lev]
      (dxh,Adxh,dxH,rH)=work_vectors[lev]
      # Pre-smooth current solution
      smooth!(smooth_caches[lev], xh, rh, Ah; maxiter=smooth_iter)

      # Restrict the residual
      mul!(rH,restrictions[lev],rh)
      # Apply next_level
      apply_GMG_level!(dxH,
                      rH,
                      lev+1,
                      mh,
                      smatrices,
                      restrictions,
                      interpolations,
                      smooth_caches,
                      work_vectors;
                      smooth_iter=smooth_iter)
      # Interpolate dxH in finer space
      mul!(dxh, interpolations[lev], dxH)
      # Update solution
      xh .= xh .+ dxh
      # Update residual
      mul!(Adxh, Ah, dxh)
      rh .= rh .- Adxh

      # Post-smooth current solution
      smooth!(smooth_caches[lev], xh, rh, Ah; maxiter=smooth_iter)

    end
  end
end


function GMG!(x::PVector,
              b::PVector,
              mh::ModelHierarchy,
              fespaces,
              smatrices,
              interpolations,
              restrictions;
              rtol=1.0e-06,
              maxiter=10,
              smooth_iter=2)

  work_vectors=allocate_work_vectors(mh,smatrices)
  smooth_caches=setup_smooth_caches(mh,smatrices)
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
                       smooth_caches,
                       work_vectors;
                       smooth_iter=smooth_iter)
    nrm_r = norm(rh)
    rel_res = nrm_r / nrm_r0
    if (GridapP4est.i_am_main(model.parts))
      @printf "%6i  %12.4e\n" current_iter rel_res
    end
    current_iter += 1
  end

  return current_iter
end
