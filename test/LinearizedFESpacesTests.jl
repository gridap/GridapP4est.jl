module LinearizedFESpacesTests
  using Gridap
  using PartitionedArrays
  using GridapDistributed
  using GridapP4est
  using FillArrays
  using Test
  using MPI

  function generate_analytical_functions(order)
    u(x)   = x[1]^order
    f(x)   = -Δ(u)(x)
    u, f  
  end

  function generate_weak_residual_form(Vh)
      function weak_res(uH, dΩ, f)
        vh_basis = Gridap.get_fe_basis(Vh)
        ∫(∇(uH)⋅∇(vh_basis))*dΩ - ∫(f*vh_basis)*dΩ
      end
  end

  function run(distribute)
    # debug_logger = ConsoleLgger(stderr, Logging.Debug)
    # global_logger(debug_logger); # Enable the debug logger globally
    ranks = distribute(LinearIndices((MPI.Comm_size(MPI.COMM_WORLD),)))

    coarse_model=CartesianDiscreteModel((0,1,0,1),(1,1))
    model=OctreeDistributedDiscreteModel(ranks,coarse_model,2)

    cell_gids=get_cell_gids(model.dmodel)
    ref_coarse_flags=map(partition(cell_gids)) do indices
      flags=Vector{Cint}(undef,local_length(indices))
      flags.=nothing_flag
      for i=1:length(indices)
         if i%2==0
            flags[i]=refine_flag
         end
      end 
      flags
    end
    modelH,_=adapt(model,ref_coarse_flags)

    for order in (2,4,8)
      degree = 2*order+1
      I(x)   = x[1]^order
      u,f=generate_analytical_functions(order)

      reffeH = ReferenceFE(lagrangian,Float64,order)
      VH     = TestFESpace(modelH,reffeH;
                          conformity=:H1,
                          dirichlet_tags="boundary")
      UH     = TrialFESpace(VH,u)

      Vh,modelh=Gridap.LinearizedFESpace(modelH,reffeH;
                                         conformity=:H1,
                                         dirichlet_tags="boundary")
      Uh=TrialFESpace(Vh,u)
      Ωh=Triangulation(modelh)
      dΩh=Measure(Ωh,degree)

      # Test L2-projection
      aL2dΩh(uH,vh)=∫(uH*vh)dΩh
      lL2Ωh(vh)    =∫(I*vh)dΩh
      op=AffineFEOperator(aL2dΩh,lL2Ωh,UH,Vh)
      Ih=solve(op)
      eh=sum(∫((I-Ih)*(I-Ih))dΩh)
      @test eh < 1.0e-12

      # Test laplacian (Petrov-Galerkin)
      adΩHh(UH,vh)=∫(∇(UH)⋅∇(vh))dΩh
      lΩHh(vh)    =∫(f*vh)dΩh
      op=AffineFEOperator(adΩHh,lΩHh,UH,Vh)
      uh=solve(op)
      eh=sum(∫((u-uh)*(u-uh))dΩh)
      @test eh < 1.0e-12

      weak_res = generate_weak_residual_form(Vh)
      rh = weak_res(uh, dΩh, f)
      rh_vec = assemble_vector(rh,Vh)
      @test norm(rh_vec) < 1.0e-12

      û(x) = sin(3.2 * x[1]^2) * cos(x[1]) + sin(4.6 * x[1]) * cos(5.2 * x[1])
      ŝ(x) = exp(x[1] / 2) + 2
      k̂(x) = 2 + sin(x[1])

      q(x) = k̂(x) * ∇(û)(x)
      f(x) = -(∇ ⋅ q)(x) + ŝ(x) * û(x)

      VH     = TestFESpace(modelH,reffeH;
                          conformity=:H1,
                          dirichlet_tags="boundary")
      UH     = TrialFESpace(VH,û)

      Vh,modelh=Gridap.LinearizedFESpace(modelH,reffeH;
                                         conformity=:H1,
                                         dirichlet_tags="boundary")

      Ωh = Triangulation(modelh)
      dΩh = Measure(Ωh, degree)
      dv = Gridap.get_fe_basis(Vh)
      a(u, v) = ∫(∇(v) ⋅ (k̂ * ∇(u)) + ŝ * u * v)dΩh
      l(v) = ∫(f * v)dΩh
      affine_op=AffineFEOperator(a, l, UH, Vh)
      ũh = solve(affine_op)

      r(u,dv)=a(u,dv)-l(dv)
      j(u,du,dv)=a(du,dv)

      tol = 1.0e-12
      op = FEOperator(r,j,UH,Vh)
      r,A = Gridap.Algebra.residual_and_jacobian(op,ũh)
      @test norm(r) < tol 
      
    end
  end
end
