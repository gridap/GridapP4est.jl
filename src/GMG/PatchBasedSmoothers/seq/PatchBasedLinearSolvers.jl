struct PatchBasedLinearSolver{A} <: Gridap.Algebra.LinearSolver
  Ph::A
end

struct PatchBasedSymbolicSetup <: Gridap.Algebra.SymbolicSetup
  solver::PatchBasedLinearSolver
end

function Gridap.Algebra.symbolic_setup(ls::PatchBasedLinearSolver,mat::AbstractMatrix)
  PatchBasedSymbolicSetup(ls)
end

mutable struct PatchBasedSmootherNumericalSetup{A,B,C,D,E} <: Gridap.Algebra.NumericalSetup
  solver         :: PatchBasedLinearSolver
  Ap             :: A
  nsAp           :: B
  rp             :: C
  dxp            :: D
  w              :: E
end

function Gridap.Algebra.numerical_setup(ss::PatchBasedSymbolicSetup,A::AbstractMatrix)
  Ph=ss.solver.Ph
  assembler=SparseMatrixAssembler(Ph,Ph)
  Ωₚ  = Triangulation(Ph.patch_decomposition)
  order=1
  dΩₚ = Measure(Ωₚ,2*order+1)
  a(u,v)=∫(∇(v)⋅∇(u))*dΩₚ
  Ap=assemble_matrix(a,assembler,Ph,Ph)
  solver = LUSolver()
  ssAp   = symbolic_setup(solver,Ap)
  nsAp   = numerical_setup(ssAp,Ap)
  rp  = zeros(size(Ap,1))
  dxp = zeros(size(Ap,2))
  w   = compute_weight_operators(ss.solver.Ph)
  PatchBasedSmootherNumericalSetup(ss.solver,Ap,nsAp,rp,dxp,w)
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
