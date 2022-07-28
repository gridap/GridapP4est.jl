using Gridap
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.ReferenceFEs
using GridapP4est
using FillArrays

# Consider a one-level smoother (multilevel comes next)
domain = (0,1,0,1)
partition = (2,2)
model = simplexify(CartesianDiscreteModel(domain,partition))

order = 2
reffe = ReferenceFE(lagrangian,Float64,order)
Vₕ    = TestFESpace(model,reffe;conformity=:H1,dirichlet_tags="boundary")
Uₕ    = TrialFESpace(Vₕ)

pd=PatchDecomposition(model)
Ph=PatchFESpace(model,reffe,H1Conformity(),pd,Vₕ)

assembler=SparseMatrixAssembler(Ph,Ph)

Ωₚ  = Triangulation(pd)
dΩₚ = Measure(Ωₚ,2*order+1)
a(u,v)=∫(∇(v)⋅∇(u))*dΩₚ
l(v)=∫(1*v)*dΩₚ
Ah=assemble_matrix(a,assembler,Ph,Ph)
fh=assemble_vector(l,assembler,Ph)
Ah\fh

x=rand(num_free_dofs(Vₕ))
y=zeros(num_free_dofs(Ph))
GridapP4est.prolongate!(y,Ph,x)
w=GridapP4est.compute_weight_operators(Ph)
x2=copy(x)
fill!(x2,0.0)
GridapP4est.inject!(x2,Ph,y)
@assert all(x .≈ x2)

Ω  = Triangulation(model)
dΩ = Measure(Ω,2*order+1)
a(u,v)=∫(∇(v)⋅∇(u))*dΩ
l(v)=∫(1.0*v)*dΩ
op=AffineFEOperator(a,l,Uₕ,Vₕ)
A=op.op.matrix
b=op.op.vector

M=PatchBasedLinearSolver(LUSolver(),Ph)

s=RichardsonSmoother(M,10,1.0/3.0)
x=zeros(num_free_dofs(Vₕ))
r=b-A*x
solve!(x,s,A,r)


# ON another note. Related to FE assembly. We are going to need:
# "Por otra parte, tb podemos tener metodos q reciben una patch-cell array y la
# aplanan para q parezca una cell array (aunq con cells repetidas). Combinando las
# patch-cell local matrices y cell_dofs aplanadas puedes usar el assembly verbatim si
# quieres ensamblar la matriz."

# Another note. During FE assembly we may end computing the cell matrix of a given cell
# more than once due to cell overlapping among patches (recall the computation of these
# matrices is lazy, it occurs on first touch). Can we live with that or should we pay
# attention on how to avoid this? I think that Gridap already includes tools for
# taking profit of this, I think it is called MemoArray, but it might be something else
# (not 100% sure, to investigate)
