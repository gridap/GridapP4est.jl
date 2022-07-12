
# By now, just simple Jacobi smoother !
# We may implement in the future block Gauss Seidel, etc.
function smooth!(x::PVector, A::PSparseMatrix, b::PVector; maxiter=1)
  r  = b - A*x   # dynamic memory alloc
  dx = copy(x)   # dynamic memory alloc
  inv_diag=map_parts(A.owned_owned_values) do a
    1.0 ./ diag(a) # dynamic memory alloc
  end
  iter=0
  while iter <= maxiter
    map_parts(inv_diag,dx.owned_values,r.owned_values) do inv_diag, dx, r
       dx .= inv_diag .* r
    end
    x .= x .+ dx
    r .= b .- A*x
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

  A      = smatrices[1]
  r      = b - A*x   # dynamic memory alloc
  nrm_r0 = norm(r)
  nrm_r  = nrm_r0
  current_iter = 0
  rel_res = nrm_r / nrm_r0
  model = get_level_model(mh,1)

  if (GridapP4est.i_am_main(model.parts))
    @printf "%6s  %12s" "Iter" "Rel res\n"
  end

  while current_iter <= maxiter && rel_res > rtol
    if (GridapP4est.i_am_main(model.parts))
       @printf "%6i  %12.4e\n" current_iter rel_res
    end
     smooth!(x, A, b; maxiter=smooth_iter)
     r .= b .- A * x  # dynamic memory alloc
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
     dx_interp = PVector(0.0,smatrices[1].cols)
     mul!(dx_interp, interpolations[1], dx)

     # Interpolate the dx
     x .= x .+ dx_interp

     smooth!(x, A, b; maxiter=smooth_iter)
     r .= b .- A * x
     nrm_r = norm(r)
     rel_res = nrm_r / nrm_r0
     current_iter += 1
  end
  return current_iter
end
