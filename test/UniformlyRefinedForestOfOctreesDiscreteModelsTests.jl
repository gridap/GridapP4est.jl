module UniformlyRefinedForestOfOctreesDiscreteModelsTests
  using MPI
  using Gridap
  using Gridap.Algebra
  using Gridap.FESpaces
  using PartitionedArrays
  using GridapDistributed
  using GridapP4est
  using Test
  using ArgParse

  function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--subdomains", "-s"
        help = "Tuple with the # of coarse subdomains per Cartesian direction"
        arg_type = Int64
        default=[1,1]
        nargs='+'
        "--num-uniform-refinements", "-r"
        help = "# of uniform refinements"
        arg_type = Int64
        default=1
    end
    return parse_args(s)
  end

  function run(parts,subdomains,num_uniform_refinements)
    # Manufactured solution
    u(x) = x[1] + x[2]
    f(x) = -Δ(u)(x)

    if length(subdomains)==2
      domain=(0,1,0,1)
    else
      @assert length(subdomains)==3
      domain=(0,1,0,1,0,1)
    end

    coarse_discrete_model=CartesianDiscreteModel(domain,subdomains)
    model=UniformlyRefinedForestOfOctreesDiscreteModel(parts,
                                                       coarse_discrete_model,
                                                       num_uniform_refinements)

    # FE Spaces
    order=1
    reffe = ReferenceFE(lagrangian,Float64,order)
    V = TestFESpace(model,reffe,dirichlet_tags="boundary")
    U = TrialFESpace(u,V)

    trian=Triangulation(model)
    dΩ=Measure(trian,2*(order+1))

    function a(u,v)
      ∫(∇(v)⋅∇(u))dΩ
    end
    function l(v)
      ∫(v*f)dΩ
    end
    dv = get_fe_basis(V)
    du = get_trial_fe_basis(U)
    assem = SparseMatrixAssembler(U,V)

    dof_values = PVector(0.0,V.gids)
    uh = FEFunction(U,dof_values)
    data = collect_cell_matrix_and_vector(U,V,a(du,dv),l(dv),uh)
    A,b = assemble_matrix_and_vector(assem,data)

    x = A\b
    r = A*x -b
    uh = FEFunction(U,x)

    # Error norms and print solution
    trian=Triangulation(model)
    dΩ=Measure(trian,2*order)
    e = u-uh
    e_l2 = sum(∫(e*e)dΩ)
    tol = 1.0e-9
    @test e_l2 < tol
    map_parts(parts) do part
     if (part==1)
       println("$(e_l2) < $(tol)\n")
     end
    end
  end
  if !MPI.Initialized()
    MPI.Init()
  end
  parsed_args = parse_commandline()
  subdomains = Tuple(parsed_args["subdomains"])
  num_uniform_refinements = parsed_args["num-uniform-refinements"]
  parts = get_part_ids(mpi,(prod(subdomains)))
  run(parts,subdomains,num_uniform_refinements)
  MPI.Finalize()
end # module
