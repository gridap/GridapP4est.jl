
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

function GMG!(x::PVector,
              b::PVector,
              mh::ModelHierarchy,
              fespaces,
              smatrices,
              mmatrices,
              interpolations,
              restrictions;
              rtol=1.0e-06,
              maxiter=10,
              smooth_iter=2)

  A = smatrices[1]
  Adx = PVector(0.0, A.cols)
  dx_interp = PVector(0.0, A.cols)
  r = PVector(0.0, A.rows)
  mul!(Adx,A,x)
  r .= b .- Adx
  nrm_r0 = norm(r)
  nrm_r  = nrm_r0
  current_iter = 0
  rel_res = nrm_r / nrm_r0
  model = get_level_model(mh,1)
  cache=smooth_cache(A)

  if (GridapP4est.i_am_main(model.parts))
    @printf "%6s  %12s" "Iter" "Rel res\n"
  end

  while current_iter <= maxiter && rel_res > rtol
    if (GridapP4est.i_am_main(model.parts))
       @printf "%6i  %12.4e\n" current_iter rel_res
    end
     # Pre-smooth current solution
     smooth!(cache, x, r, A; maxiter=smooth_iter)

     nrm_r = norm(r)
     rel_res = nrm_r / nrm_r0
     if (GridapP4est.i_am_main(model.parts))
       @printf "%6i  %12.4e\n" current_iter rel_res
     end
     modelH = get_level_model(mh,2)
     if (GridapP4est.i_am_in(modelH.parts))
       AH = smatrices[2]
       y  = PVector(0.0,AH.cols)
     else
       y=nothing
     end

     # Restrict the residual
     mul!(y,restrictions[1],r)

     if (GridapP4est.i_am_in(modelH.parts))
      dx=PVector(0.0,smatrices[2].cols)
      # TO-DO: replace this coarsest solver by a call to a
      #        a sparse direct solver, e.g., PETSc+MUMPs
      map_parts(smatrices[2].owned_owned_values,dx.owned_values,y.owned_values) do A, dx, y
        dx .= A \ y
      end
      #IterativeSolvers.cg!(dx,
      #                     smatrices[2],
      #                     y;
      #                     verbose=true,
      #                     reltol=1.0e-14)
     else
      dx=nothing
     end
     # Interpolate the dx
     mul!(dx_interp, interpolations[1], dx)

     # Update solution
     x .= x .+ dx_interp

     # Update residual
     mul!(Adx, A, dx_interp)
     r .= r .- Adx

     # Post-smooth current solution
     smooth!(cache, x, r, A; maxiter=smooth_iter)

     nrm_r = norm(r)
     rel_res = nrm_r / nrm_r0
     current_iter += 1
  end
  return current_iter
end
