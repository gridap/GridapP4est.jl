
struct JacobiLinearSolver <: Gridap.Algebra.LinearSolver
end

struct JacobiSymbolicSetup <: Gridap.Algebra.SymbolicSetup
end

function Gridap.Algebra.symbolic_setup(jacobi_smoother::JacobiLinearSolver,mat::AbstractMatrix)
  JacobiSymbolicSetup()
end

mutable struct JacobiNumericalSetup{A} <: Gridap.Algebra.NumericalSetup
  inv_diag :: A
end

function Gridap.Algebra.numerical_setup(ss::JacobiSymbolicSetup,a::AbstractMatrix)
  inv_diag=1.0./diag(a)
  JacobiNumericalSetup(inv_diag)
end

function Gridap.Algebra.numerical_setup!(ns::JacobiNumericalSetup, A::AbstractMatrix)
  ns.inv_diag .= 1.0 ./ diag(a)
end

function Gridap.Algebra.numerical_setup(ss::JacobiSymbolicSetup,A::PSparseMatrix)
  inv_diag=map_parts(A.owned_owned_values) do a
    1.0 ./ diag(a)
  end
  JacobiNumericalSetup(inv_diag)
end

function Gridap.Algebra.numerical_setup!(ns::JacobiNumericalSetup, A::PSparseMatrix)
  map_parts(ns.inv_diag,A.owned_owned_values) do inv_diag, a
    inv_diag .= 1.0 ./ diag(a)
  end
  ns
end

function Gridap.Algebra.solve!(
  x::AbstractVector,ns::JacobiNumericalSetup,r::AbstractVector)
  inv_diag=ns.inv_diag
  x .= inv_diag .* r
  x
end

function Gridap.Algebra.solve!(
  x::PVector,ns::JacobiNumericalSetup,r::PVector)
  inv_diag=ns.inv_diag
  map_parts(inv_diag,x.owned_values,r.owned_values) do inv_diag, x, r
    x .= inv_diag .* r
  end
  x
end
