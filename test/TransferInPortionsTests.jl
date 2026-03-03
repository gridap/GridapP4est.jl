module TransferInPortionsTests
  using P4est_wrapper
  using GridapP4est
  using Gridap
  using PartitionedArrays
  using GridapDistributed
  using MPI
  using Gridap.FESpaces
  using Gridap.CellData
  using FillArrays
  using Logging

  include("CoarseDiscreteModelsTools.jl")

  function _generate_triangulation_portion(ranks, fmodel)
    trians = map(ranks, 
         local_views(fmodel.dmodel), 
                     partition(get_cell_gids(fmodel))) do rank, lmodel, indices
        mask = Vector{Bool}(undef,num_cells(lmodel))
        mask .= false
        cell_to_part = local_to_owner(indices)
        graph = GridapDistributed.compute_cell_graph(lmodel)
        ncells = num_cells(lmodel)
        icell_to_jcells_ptrs = graph.colptr
        icell_to_jcells_data = graph.rowval
        for icell in 1:ncells
           if cell_to_part[icell] == rank
              pini = icell_to_jcells_ptrs[icell]
              pend = icell_to_jcells_ptrs[icell+1]-1
              for p in pini:pend
                 jcell = icell_to_jcells_data[p]
                 if cell_to_part[jcell] != rank
                    mask[icell] = true
                 end                 
              end
           end
        end
        Triangulation(lmodel, mask)
    end
    GridapDistributed.DistributedTriangulation(trians,fmodel)
  end

  function _generate_triangulation_portion(ranks, fmodel, ctrian, glue)
    trians = map(ranks,
             local_views(fmodel.dmodel),
             local_views(ctrian),
             glue) do rank, fmodel, ctrian, glue
        mask = Vector{Bool}(undef,num_cells(fmodel))
        mask .= false

        cglue = get_glue(ctrian, Val{num_cell_dims(ctrian)}())
        for ccell in cglue.tface_to_mface
            fcells = glue.o2n_faces_map[ccell]
            for fcell in fcells
                mask[fcell] = true
            end
        end
        Triangulation(fmodel, mask)
    end
    GridapDistributed.DistributedTriangulation(trians,fmodel)
  end

  function generate_triangulation_portion(ranks,fmodel; ctrian=nothing, glue=nothing)
    if (length(ranks)==1 && ctrian==nothing)
       trians = map(ranks, local_views(fmodel.dmodel)) do rank, fmodel
            mask = Vector{Bool}(undef,num_cells(fmodel))
            mask .= false
            for cell=1:Int(round(num_cells(fmodel)*0.25))
               mask[cell] = true
            end
            for cell=Int(round(num_cells(fmodel)*0.75)):num_cells(fmodel)
               mask[cell] = true
            end
            Triangulation(fmodel, mask)
        end
        return GridapDistributed.DistributedTriangulation(trians,fmodel)
    else 
        if ctrian==nothing
            @assert glue==nothing
            return _generate_triangulation_portion(ranks, fmodel)
        else
            return _generate_triangulation_portion(ranks, fmodel, ctrian, glue)
        end
    end
  end

  function _generate_coarsened_triangulation_portion(ranks, gmodel)
    trians = map(local_views(gmodel.dmodel)) do gmodel
        mask = Vector{Bool}(undef,num_cells(gmodel))
        mask .= false
        # Select first quarter and last quarter of cells
        for cell=1:Int(round(num_cells(gmodel)*0.25))
           mask[cell] = true
        end
        for cell=Int(round(num_cells(gmodel)*0.75)):num_cells(gmodel)
           mask[cell] = true
        end
        Triangulation(gmodel, mask)
    end
    return GridapDistributed.DistributedTriangulation(trians,gmodel)
  end

  function generate_analytical_problem_functions(T::Type{Float64},order)
    u(x) = x[1]+x[2]^order
    f(x) = -Δ(u)(x)
    u,f
  end

  function generate_analytical_problem_functions(T::Type{VectorValue{2,Float64}},order)
    u(x) = VectorValue(x[1]+x[2]^order,x[1]^order+x[2])
    f(x) = -Δ(u)(x)
    u,f
  end

  function test_fe_space_on_triangulation(ranks,cmodel,ctrian,order,T::Type,amr_step)
    
    u,_ = generate_analytical_problem_functions(T,order)
 
    # FINE TRIANGULATION (1)
    flags = map(ranks,partition(get_cell_gids(cmodel.dmodel))) do rank,indices
        flags=zeros(Cint,length(indices))
        flags.=nothing_flag
        
        flags[1]=refine_flag
        flags[own_length(indices)]=refine_flag
        
        # To create some unbalance
        if (rank%2==0 && own_length(indices)>1)
              flags[div(own_length(indices),2)]=refine_flag
        end
        flags
    end
    fmodel,glue=Gridap.Adaptivity.adapt(cmodel,flags)
    ftrian = generate_triangulation_portion(ranks,fmodel,ctrian=ctrian,glue=glue)

    reffe = ReferenceFE(lagrangian,T,order)
    Vh1 = FESpace(ftrian,reffe,conformity=:H1)
    Uh1 = TrialFESpace(Vh1,u)
    uh1 = interpolate(u,Uh1)
    
    # COARSENED TRIANGULATION (2)
    flags = map(ranks,partition(get_cell_gids(fmodel))) do rank,indices
        flags = zeros(Cint,length(indices))
        flags .= nothing_flag
        n = own_length(indices)
        if n >= 4
          flags[max(1,n-3):n] .= coarsen_flag
        end
        flags
    end
    gmodel,_ = Gridap.Adaptivity.adapt(fmodel,flags)
    gtrian = _generate_coarsened_triangulation_portion(ranks, gmodel)

    Vh2 = FESpace(gtrian,reffe,conformity=:H1)
    Uh2 = TrialFESpace(Vh2,u)
    uh2 = interpolate(u,Uh2)

    # Interpolate fine FE function on the coarsened triangulation
    uh1_on_Uh2 = interpolate(uh1,Uh2) 

    # Check that the interpolation is accurate
    e = u - uh1_on_Uh2
    Ωg = gtrian
    dΩg = Measure(Ωg, 2*order+1)
    el2 = sqrt(sum(∫(e⋅e)*dΩg))
    tol = 1.0e-10
    @assert el2 < tol "Transfer error too large: el2 = $el2"

    fmodel, ftrian
  end
 
  function test_2d_fe_space_on_triangulation(ranks,order,T::Type;num_amr_steps=3,num_ghost_layers=1)
    coarse_model=CartesianDiscreteModel((0,1,0,1),(1,1))
    dmodel=OctreeDistributedDiscreteModel(ranks,coarse_model,2;num_ghost_layers=num_ghost_layers)
    dtrian=generate_triangulation_portion(ranks,dmodel)
    for i=1:num_amr_steps
     dmodel,dtrian=test_fe_space_on_triangulation(ranks,dmodel,dtrian,order,T,i)
    end
  end 
  
  function _field_type(::Val{Dc}, scalar_or_vector::Symbol) where Dc
    if scalar_or_vector==:scalar
      Float64
    else 
      @assert scalar_or_vector==:vector
      VectorValue{Dc,Float64}
    end
  end 

  function run(distribute)
    ranks = distribute(LinearIndices((MPI.Comm_size(MPI.COMM_WORLD),)))
    for order=1:2, scalar_or_vector in (:scalar,), num_ghost_layers in (1,)
       test_2d_fe_space_on_triangulation(ranks,
                                         order,
                                         _field_type(Val{2}(),scalar_or_vector),
                                         num_amr_steps=3,
                                         num_ghost_layers=num_ghost_layers)
    end
  end
end
