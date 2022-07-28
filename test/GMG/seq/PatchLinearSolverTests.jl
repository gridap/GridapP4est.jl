using Gridap
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.ReferenceFEs
using GridapP4est
using FillArrays
using PartitionedArrays
using Test
const order=2

function returns_Ph_xh_Vh(model)
  reffe = ReferenceFE(lagrangian,Float64,order)
  Vh    = TestFESpace(model,reffe;conformity=:H1)
  PD=PatchDecomposition(model)
  Ph=PatchFESpace(model,reffe,H1Conformity(),PD,Vh)
  assembler=SparseMatrixAssembler(Ph,Ph)
  Ωₚ  = Triangulation(PD)
  dΩₚ = Measure(Ωₚ,2*order+1)
  a(u,v)=∫(∇(v)⋅∇(u))*dΩₚ
  l(v)=∫(1*v)*dΩₚ
  Ah=assemble_matrix(a,assembler,Ph,Ph)
  fh=assemble_vector(l,assembler,Ph)
  Ph,Ah\fh,Vh
end

domain=(0.0,1.0,0.0,1.0)
partition=(2,3)
parts=get_part_ids(sequential,(1,2))
dmodel=CartesianDiscreteModel(parts,domain,partition)
model=CartesianDiscreteModel(domain,partition)
dPh,dxh,dVh=returns_Ph_xh_Vh(dmodel);
Ph,xh,Vh=returns_Ph_xh_Vh(model)
@test num_free_dofs(Ph) == num_free_dofs(dPh)
@test all(dxh.owned_values.parts[1] .≈ xh[1:5])
@test all(dxh.owned_values.parts[2] .≈ xh[6:end])

function compute_matrix_vector(model,Vh)
  Ω = Triangulation(model)
  dΩ = Measure(Ω,2*order+1)
  a(u,v)=∫(∇(v)⋅∇(u))*dΩ
  l(v)=∫(1*v)*dΩ
  assembler=SparseMatrixAssembler(Vh,Vh)
  Ah=assemble_matrix(a,assembler,Vh,Vh)
  lh=assemble_vector(l,assembler,Vh)
  Ah,lh
end



function test_smoother(Ph,Vh,A,b)
  M=PatchBasedLinearSolver(LUSolver(),Ph)
  s=RichardsonSmoother(M,10,1.0/3.0)
  x=GridapP4est._allocate_col_vector(A)
  println(typeof(x))
  r=b-A*x
  solve!(x,s,A,r)
  x
end

domain=(0.0,1.0,0.0,1.0)
partition=(2,3)
parts=get_part_ids(sequential,(1,1))
dmodel=CartesianDiscreteModel(parts,domain,partition)
model=CartesianDiscreteModel(domain,partition)
dPh,dxh,dVh=returns_Ph_xh_Vh(dmodel);
Ph,xh,Vh=returns_Ph_xh_Vh(model)
A,b=compute_matrix_vector(model,Vh)
dA,db=compute_matrix_vector(dmodel,dVh);

x=test_smoother(Ph,Vh,A,b)
dx=test_smoother(dPh,dVh,dA,db)
