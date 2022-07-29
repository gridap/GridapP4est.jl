struct PatchBasedLinearSolver{A,B} <: Gridap.Algebra.LinearSolver
  PD :: A
  Ph :: B
  M  :: Gridap.Algebra.LinearSolver
end

struct PatchBasedSymbolicSetup <: Gridap.Algebra.SymbolicSetup
  solver :: PatchBasedLinearSolver
end

function Gridap.Algebra.symbolic_setup(ls::PatchBasedLinearSolver,mat::AbstractMatrix)
  PatchBasedSymbolicSetup(ls)
end

struct PatchBasedSmootherNumericalSetup{A,B,C,D,E} <: Gridap.Algebra.NumericalSetup
  solver         :: PatchBasedLinearSolver
  Ap             :: A
  nsAp           :: B
  rp             :: C
  dxp            :: D
  w              :: E
end

function Gridap.Algebra.numerical_setup(ss::PatchBasedSymbolicSetup,A::AbstractMatrix)
  PD=ss.solver.PD
  Ph=ss.solver.Ph
  assembler=SparseMatrixAssembler(Ph,Ph)
  Ωₚ  = Triangulation(PD)
  order=1
  dΩₚ = Measure(Ωₚ,2*order+1)
  a(u,v)=∫(∇(v)⋅∇(u))*dΩₚ
  Ap=assemble_matrix(a,assembler,Ph,Ph)
  solver = ss.solver.M
  ssAp   = symbolic_setup(solver,Ap)
  nsAp   = numerical_setup(ssAp,Ap)
  rp  = _allocate_row_vector(Ap)
  dxp = _allocate_col_vector(Ap)
  w   = compute_weight_operators(Ph)
  PatchBasedSmootherNumericalSetup(ss.solver,Ap,nsAp,rp,dxp,w)
end

function _allocate_col_vector(A::AbstractMatrix)
  zeros(size(A,2))
end

function _allocate_row_vector(A::AbstractMatrix)
  zeros(size(A,1))
end

function _allocate_col_vector(A::PSparseMatrix)
  PVector(0.0,A.cols)
end

function _allocate_row_vector(A::PSparseMatrix)
  PVector(0.0,A.rows)
end

function Gridap.Algebra.numerical_setup!(ns::PatchBasedSmootherNumericalSetup, A::AbstractMatrix)
  Gridap.Helpers.@notimplemented
end

function Gridap.Algebra.solve!(
  x::AbstractVector,ns::PatchBasedSmootherNumericalSetup,r::AbstractVector)
  Ap,nsAp,rp,dxp,w=ns.Ap,ns.nsAp,ns.rp,ns.dxp,ns.w
  prolongate!(rp,ns.solver.Ph,r)
  solve!(dxp,nsAp,rp)
  inject!(x,ns.solver.Ph,dxp,w)
end
