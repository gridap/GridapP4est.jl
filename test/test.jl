using Gridap
using PartitionedArrays
using GridapDistributed
using GridapP4est
using MPI
using FillArrays
using Test

include("CoarseDiscreteModelsTools.jl")

MPI.Init()
nprocs = MPI.Comm_size(MPI.COMM_WORLD)
ranks  = with_mpi() do distribute
  distribute(LinearIndices((prod(nprocs),)))
end

coarse_model=setup_model(Val{2},1)
dmodel=OctreeDistributedDiscreteModel(ranks,coarse_model,1)

ref_coarse_flags=map(ranks,partition(get_cell_gids(dmodel.dmodel))) do rank,indices
    flags=zeros(Int,length(indices))
    flags.=nothing_flag        
    flags[2]=refine_flag
    flags
end

fmodel,glue=adapt(dmodel,ref_coarse_flags);

Λ=SkeletonTriangulation(fmodel)
dΛ=Measure(Λ,2)

Q = FESpace(fmodel,
            ReferenceFE(lagrangian,Float64,1); 
            conformity=:L2)

qh_dofs=pfill(1.0, partition(Q.gids))
qh=FEFunction(Q,qh_dofs)

dcm=∫(qh.minus)dΛ
dcp=∫(qh.plus)dΛ


Ω = Triangulation(fmodel)
map(local_views(∫(qh.minus)dΛ),local_views(Λ),local_views(Ω)) do dc,Λ,Ω
  println(get_array(dc)) 
  Gridap.CellData.is_change_possible(Λ,Ω)
end 

 ∫(qh.plus)dΛ

Λd=SkeletonTriangulation(dmodel)
dΛd=Measure(Λd,2)

Qd = FESpace(dmodel,
            ReferenceFE(lagrangian,Float64,1); 
            conformity=:L2)

qhd_dofs=pfill(1.0, partition(Qd.gids))
qhd=FEFunction(Qd,qhd_dofs)

∫(qhd.plus)dΛd

# Now with DG
h = 2
γ = 10


Ω = Triangulation(fmodel)
#Γn = Boundary(fmodel,tags="neumann")
#n_Γn = get_normal_vector(Γn)

k = 2
dΩ = Measure(Ω,2*k)
u((x,y)) = (x+y)^k
f(x) = -Δ(u,x)
#g = n_Γn⋅∇(u)

reffe = ReferenceFE(lagrangian,Float64,k)
V_dg = FESpace(fmodel,reffe,conformity=:L2)

Λ = Skeleton(fmodel)
Γd = Boundary(fmodel,tags="boundary")

dΛ = Measure(Λ,2*k)
dΓd = Measure(Γd,2*k)

n_Γd = get_normal_vector(Γd)
n_Λ = get_normal_vector(Λ)

v=get_fe_basis(V_dg)
ub=get_trial_fe_basis(V_dg)

ub*n_Λ

(γ/h)*jump(v*n_Λ)⋅jump(ub*n_Λ)

a_dg(u,v) =
  ∫( ∇(v)⋅∇(u) )*dΩ +
  ∫( (γ/h)*v*u  - v*(n_Γd⋅∇(u)) - (n_Γd⋅∇(v))*u )*dΓd +
  ∫( (γ/h)*jump(v*n_Λ)⋅jump(u*n_Λ) -
     jump(v*n_Λ)⋅mean(∇(u)) -
     mean(∇(v))⋅jump(u*n_Λ) )*dΛ

l_dg(v) =
  ∫( v*f )*dΩ +
  #∫( v*g )dΓn +
  ∫( (γ/h)*v*u - (n_Γd⋅∇(v))*u )*dΓd

op = AffineFEOperator(a_dg,l_dg,V_dg,V_dg)
uh = solve(op)
eh = u - uh
@test sqrt(sum( ∫(abs2(eh))dΩ )) < 1.0e-9