
struct JacobiSmoother <: Gridap.Algebra.LinearSolver
  num_smooth_steps::Int
end

struct JacobiSmootherSymbolicSetup <: Gridap.Algebra.SymbolicSetup
  jacobi_smoother::JacobiSmoother
end

function Gridap.Algebra.symbolic_setup(jacobi_smoother::JacobiSmoother,mat::PSparseMatrix)
  JacobiSmootherSymbolicSetup(jacobi_smoother)
end

mutable struct JacobiSmootherNumericalSetup{A,B,C,D} <: Gridap.Algebra.NumericalSetup
  jacobi_smoother:: JacobiSmoother
  A              :: A
  Adx            :: B
  dx             :: C
  inv_diag       :: D
end

function Gridap.Algebra.numerical_setup(ss::JacobiSmootherSymbolicSetup,A::PSparseMatrix)
  Adx = PVector(0.0,A.rows)
  dx  = PVector(0.0,A.cols)
  inv_diag=map_parts(A.owned_owned_values) do a
    1.0 ./ diag(a)
  end
  JacobiSmootherNumericalSetup(ss.jacobi_smoother,A,Adx,dx,inv_diag)
end

function Gridap.Algebra.numerical_setup!(ns::JacobiSmootherNumericalSetup, A::PSparseMatrix)
  map_parts(ns.inv_diag,A.owned_owned_values) do inv_diag, a
    inv_diag .= 1.0 ./ diag(a)
  end
  ns.A=A
  ns
end

function Gridap.Algebra.solve!(
  x::PVector,ns::JacobiSmootherNumericalSetup,r::PVector)

  Adx,dx,inv_diag=ns.Adx,ns.dx,ns.inv_diag

  iter=0
  while iter <= ns.jacobi_smoother.num_smooth_steps
    map_parts(inv_diag,dx.owned_values,r.owned_values) do inv_diag, dx, r
       dx .= inv_diag .* r
    end
    x .= x .+ dx
    mul!(Adx, ns.A, dx)
    r .= r .- Adx
    iter += 1
  end
end
