struct PatchBasedSmoother <: Gridap.Algebra.LinearSolver
  num_smooth_steps::Int
  Ph::PatchFESpace
end

struct PatchBasedSymbolicSetup <: Gridap.Algebra.SymbolicSetup
  smoother::PatchBasedSmoother
end

function Gridap.Algebra.symbolic_setup(smoother::PatchBasedSmoother,mat::AbstractMatrix)
  PatchBasedSymbolicSetup(smoother)
end

mutable struct PatchBasedSmootherNumericalSetup{A,B,C,D,E,F,G,H} <: Gridap.Algebra.NumericalSetup
  smoother       :: PatchBasedSmoother
  A              :: A
  Adx            :: B
  dx             :: C
  Ap             :: D
  nsAp           :: E
  rp             :: F
  dxp            :: G
  w              :: H
end

function Gridap.Algebra.numerical_setup(ss::PatchBasedSymbolicSetup,A::AbstractMatrix)
  Adx = zeros(size(A,1))  # PVector(0.0,A.rows)
  dx  = zeros(size(A,2))  # PVector(0.0,A.cols)

  Ph=ss.smoother.Ph
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
  w   = compute_weight_operators(ss.smoother.Ph)
  PatchBasedSmootherNumericalSetup(ss.smoother,A,Adx,dx,Ap,nsAp,rp,dxp,w)
end

function Gridap.Algebra.numerical_setup!(ns::PatchBasedSmootherNumericalSetup, A::AbstractMatrix)
  Gridap.Helpers.@notimplemented
end

function Gridap.Algebra.solve!(
  x::AbstractVector,ns::PatchBasedSmootherNumericalSetup,r::AbstractVector)

  Adx,dx,Ap,nsAp,rp,dxp,w=ns.Adx,ns.dx,ns.Ap,ns.nsAp,ns.rp,ns.dxp,ns.w

  iter=0
  while iter <= ns.smoother.num_smooth_steps
    prolongate!(rp,ns.smoother.Ph,r)
    solve!(dxp,nsAp,rp)
    inject!(dx,ns.smoother.Ph,dxp,w)
    x .= x .+ dx
    mul!(Adx, ns.A, dx)
    r .= r .- Adx
    println(norm(r))
    iter += 1
  end
end
