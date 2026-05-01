# Run as: mpirun -np 2 julia --project=. mwe.jl
module MWE

using P4est_wrapper
using GridapP4est
using Gridap
using PartitionedArrays
using GridapDistributed
using MPI
using Gridap.FESpaces
using FillArrays
using Logging
using Test
using MPI
using PartitionedArrays

u_ex_2D((x,y)) = 2*VectorValue(-y,x)
f_ex_2D(x) = u_ex_2D(x)
u_ex_3D((x,y,z)) = 2*VectorValue(-y,x,0.) - VectorValue(0.,-z,y)
f_ex_3D(x) = u_ex_3D(x)

function get_analytical_functions(Dc)
  if Dc==2
    return u_ex_2D, f_ex_2D
  else
    @assert Dc==3
    return u_ex_3D, f_ex_3D
  end
end

function solve_maxwell(model::GridapDistributed.DistributedDiscreteModel{Dc},order) where {Dc}
  u_ex, f_ex=get_analytical_functions(Dc)

  reffe  = ReferenceFE(nedelec,order)

  V = FESpace(model.dmodel,
              reffe,
              conformity=:Hcurl,
              dirichlet_tags="boundary")


  # cell_reffe = map(local_views(model)) do model
  #    ReferenceFE(model,reffe)
  # end
  # sign_flip = Gridap.FESpaces.get_sign_flip(model, cell_reffe)
  # map(sign_flip) do sign_flip
  #   if (MPI.Comm_rank(MPI.COMM_WORLD)==0)  
  #      println("sign flip:"); println(sign_flip[1])
  #   end
  # end
  
  U = TrialFESpace(V,u_ex)
  
  trian = Triangulation(model)
  degree = 2*(order+1)
  dΩ = Measure(trian,degree)
      
  a(u,v) = ∫( (∇×u)⋅(∇×v) + u⋅v)dΩ

  dv = get_fe_basis(V)
  du = get_trial_fe_basis(U)
  dc = a(du,dv)

  map(local_views(∇×(dv))) do curl_dv
      if MPI.Comm_rank(MPI.COMM_WORLD)==1
            curl_dv_array = Gridap.CellData.get_data(curl_dv)
            print_op_tree(curl_dv_array)
            println("curl dv:"); println(evaluate(curl_dv_array[1],[Point(0.5,0.5,0.5)]))
      end
  end

  # map(local_views(V)) do local_model
  #     if MPI.Comm_rank(MPI.COMM_WORLD)==0
  #          println("local model:"); println(get_cell_dof_ids(local_model))
  #     end
  # end

  # map(local_views(dc)) do dc
  #     if MPI.Comm_rank(MPI.COMM_WORLD)==1
  #          println("local stiffness matrix:"); println(get_array(dc))
  #     end
  # end

  l(v) = ∫(f_ex⋅v)dΩ

  op = AffineFEOperator(a,l,U,V)
  if (num_free_dofs(U)==0)
    # UMFPACK cannot handle empty linear systems
    uh = zero(U)
  else
    uh = solve(op)
  end
  # uh=interpolate(u_ex,U)
  
  # uh_ex=interpolate(u_ex_3D,U)
  # map(local_views(get_free_dof_values(uh_ex)), local_views(op.op.matrix), local_views(op.op.vector)) do U_ex, A, b 
  #   r_ex = A*U_ex - b
  #   println(r_ex)
  # end
  uh,U
end 

function check_error_maxwell(model::GridapDistributed.DistributedDiscreteModel{Dc},order,uh) where {Dc}
  trian = Triangulation(model)
  degree = 2*(order+1)
  dΩ = Measure(trian,degree)

  u_ex, f_ex = get_analytical_functions(Dc)
  
  # map(local_views(uh)) do uh
  #     if MPI.Comm_rank(MPI.COMM_WORLD)==0
  #          print(get_cell_dof_values(uh)); print("\n")
  #     end
  # end
  eu = u_ex - uh

  l2(v) = sqrt(sum(∫(v⋅v)*dΩ))
  hcurl(v) = sqrt(sum(∫(v⋅v + (∇×v)⋅(∇×v))*dΩ))
  
  eu_l2 = l2(eu)
  eu_hcurl = hcurl(eu)
  
  tol = 1.0e-6
  @test eu_l2 < tol
  @test eu_hcurl < tol
end 

function run(distribute)
  np    = MPI.Comm_size(MPI.COMM_WORLD)
  ranks = distribute(LinearIndices((np,)))    
  coarse_model=CartesianDiscreteModel((0,1,0,1,0,1),(2,1,1))
  dmodel=OctreeDistributedDiscreteModel(ranks,coarse_model,0)
  order=1
  uH,UH=solve_maxwell(dmodel,order)
  check_error_maxwell(dmodel,order,uH)
end

if !MPI.Initialized()
  MPI.Init()
end

with_mpi() do distribute 
  run(distribute)
end 

MPI.Finalize()
end # module