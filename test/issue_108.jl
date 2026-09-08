
using MPI
using PartitionedArrays

using Gridap
using Gridap.Adaptivity
using Gridap.Arrays
using Gridap.Geometry

using GridapDistributed
using GridapDistributed: DistributedDiscreteModel, DistributedTriangulation

using GridapP4est
using GridapP4est: OctreeDistributedDiscreteModel, refine_flag, nothing_flag
using P4est_wrapper

Gridap.Geometry.Triangulation(t::Gridap.Adaptivity.AdaptedTriangulation,
                              mask::AbstractArray{<:Bool}) =
    Gridap.Geometry.Triangulation(t.trian, mask)

function mwe_sphere_3d_hanging_glue_on_portion(distribute, nparts::Int, n::Int)
    ranks = distribute(LinearIndices((nparts,)))
    rank  = MPI.Comm_rank(MPI.COMM_WORLD)

    GridapP4est.with(ranks; p4est_verbosity_level=P4est_wrapper.SC_LP_ERROR) do

    dim   = 3
    order = 2

    domain = (0.5, 1.0, 0.5, 1.0, 0.5, 1.0)
    cells  = (n, n, n)

    outdir = joinpath(@__DIR__, "mwe_sphere_out_no_embedded")
    rank == 0 && mkpath(outdir)
    MPI.Barrier(MPI.COMM_WORLD)

    coarse_model = CartesianDiscreteModel(domain, cells)
    model = OctreeDistributedDiscreteModel(ranks, coarse_model, 0)

    flags = map(partition(get_cell_gids(model)), ranks) do indices, r
        f = fill(Int(nothing_flag), local_length(indices))
        if r == 2
            f[2] = Int(refine_flag)
        end
        f
    end
    fmodel, _ = Gridap.Adaptivity.adapt(model, flags)
    Ωh = Triangulation(fmodel)

    # writevtk(Ωh, joinpath(outdir, "adapted_mesh"))

    cell_active = map(partition(get_cell_gids(fmodel)), ranks) do indices, r
        f = falses(local_length(indices))
        o2l = own_to_local(indices)
        if r == 1
            f[o2l[[1,2,3,4,5,7]]] .= true
        elseif r == 2
            f[o2l[[1,2,3,4,5,6,8,11,12]]] .= true
        elseif r == 3
            f[o2l[[1]]] .= true
        end
        f
    end

    Ω1_trians = map(local_views(Ωh), cell_active) do trian, mask
        Triangulation(trian, mask)
    end
    Ω1 = DistributedTriangulation(Ω1_trians, fmodel)
    dΩ1 = Measure(Ω1, 2*order)

    reffe_u = ReferenceFE(lagrangian, VectorValue{dim,Float64}, order+1)

    # writevtk(Ω1, joinpath(outdir, "Ω1"))
    # println("  V1 on Ω1 ..."); flush(stdout)

    V1 = TestFESpace(Ω1, reffe_u, conformity=:H1) # test for nparts = 3

    u(x) = VectorValue(x[3]^2,x[2]^2+x[3]^2,x[1]^2+x[2]^2)
    uh = interpolate_everywhere(u, V1) # At ab984f3e raises the error:

# ┌ Error: 
# │   exception =
# │    ArgumentError: column indices J[k] must satisfy 1 <= J[k] <= n
# │    Stacktrace:
# │      [1] sparse!(I::Vector{Int32}, J::Vector{Int32}, V::Vector{Int32}, m::Int64, n::Int64, combine::typeof(+), klasttouch::Vector{Int32}, csrrowptr::Vector{Int32}, csrcolval::Vector{Int32}, csrnzval::Vector{Int32}, csccolptr::Vector{Int32}, cscrowval::Vector{Int32}, cscnzval::Vector{Int32})
# │        @ SparseArrays ~/.julia/juliaup/julia-1.12.6+0.x64.linux.gnu/share/julia/stdlib/v1.12/SparseArrays/src/sparsematrix.jl:1205
# │      [2] sparse(I::Vector{Int32}, J::Vector{Int32}, V::Vector{Int32}, m::Int64, n::Int64, combine::Function)
# │        @ SparseArrays ~/.julia/juliaup/julia-1.12.6+0.x64.linux.gnu/share/julia/stdlib/v1.12/SparseArrays/src/sparsematrix.jl:1105
# │      [3] sparse
# │        @ ~/.julia/juliaup/julia-1.12.6+0.x64.linux.gnu/share/julia/stdlib/v1.12/SparseArrays/src/sparsematrix.jl:1340 [inlined]
# │      [4] (::PartitionedArrays.var"#find_rcv_ids_gather_scatter##0#find_rcv_ids_gather_scatter##1")(snd_ids_main::JaggedArray{Int32, Int32})
# │        @ PartitionedArrays ~/.julia/packages/PartitionedArrays/MVmxR/src/primitives.jl:651
# │      [5] map
# │        @ ~/.julia/packages/PartitionedArrays/MVmxR/src/mpi_array.jl:221 [inlined]
# │      [6] find_rcv_ids_gather_scatter(snd_ids::MPIArray{Vector{Int32}, 1})
# │        @ PartitionedArrays ~/.julia/packages/PartitionedArrays/MVmxR/src/primitives.jl:638
# │      [7] ExchangeGraph_impl_with_find_rcv_ids
# │        @ ~/.julia/packages/PartitionedArrays/MVmxR/src/primitives.jl:630 [inlined]
# │      [8] #ExchangeGraph#31
# │        @ ~/.julia/packages/PartitionedArrays/MVmxR/src/primitives.jl:594 [inlined]
# │      [9] ExchangeGraph
# │        @ ~/.julia/packages/PartitionedArrays/MVmxR/src/primitives.jl:585 [inlined]
# │     [10] compute_assembly_neighbors(indices::MPIArray{LocalIndices, 1}; kwargs::@Kwargs{})
# │        @ PartitionedArrays ~/.julia/packages/PartitionedArrays/MVmxR/src/p_range.jl:411
# │     [11] compute_assembly_neighbors
# │        @ ~/.julia/packages/PartitionedArrays/MVmxR/src/p_range.jl:399 [inlined]
# │     [12] assembly_neighbors(indices::MPIArray{LocalIndices, 1}; kwargs::@Kwargs{})
# │        @ PartitionedArrays ~/.julia/packages/PartitionedArrays/MVmxR/src/p_range.jl:386
# │     [13] assembly_neighbors
# │        @ ~/.julia/packages/PartitionedArrays/MVmxR/src/p_range.jl:380 [inlined]
# │     [14] p_vector_cache_impl(::Type, vector_partition::MPIArray{Vector{Float64}, 1}, index_partition::MPIArray{LocalIndices, 1})
# │        @ PartitionedArrays ~/.julia/packages/PartitionedArrays/MVmxR/src/p_vector.jl:291
# │     [15] p_vector_cache
# │        @ ~/.julia/packages/PartitionedArrays/MVmxR/src/p_vector.jl:259 [inlined]
# │     [16] PVector
# │        @ ~/.julia/packages/PartitionedArrays/MVmxR/src/p_vector.jl:185 [inlined]
# │     [17] PVector
# │        @ ~/.julia/packages/PartitionedArrays/MVmxR/src/p_vector.jl:523 [inlined]
# │     [18] _find_vector_type(spaces::MPIArray{Gridap.FESpaces.FESpaceWithLinearConstraints{Gridap.FESpaces.UnconstrainedFESpace{Vector{Float64}, Gridap.FESpaces.CompressedCellConformity{Vector{Int8}}}}, 1}, gids::PRange{MPIArray{LocalIndices, 1}}; split_own_and_ghost::Bool)
# │        @ GridapDistributed ~/.julia/packages/GridapDistributed/O0Pts/src/FESpaces.jl:624
# │     [19] _add_constraints(pXest_refinement_rule_type::GridapP4est.PXestUniformRefinementRuleType, models::MPIArray{UnstructuredDiscreteModel{3, 3, Float64, NonOriented}, 1}, non_conforming_glue::MPIArray{GridapP4est.NonConformingGlue{3, Vector{Int64}, Vector{Int64}, Vector{Vector{Tuple{Int64, Int64, Int64}}}, Vector{Vector{Int64}}, Vector{Vector{Int64}}, Vector{Vector{Int64}}, Vector{Dict{Int64, Tuple{Int64, Int64, Int64}}}}, 1}, trian::DistributedTriangulation{3, 3, MPIArray{AdaptedTriangulation{3, 3, BodyFittedTriangulation{3, 3, UnstructuredDiscreteModel{3, 3, Float64, NonOriented}, GridPortion{3, 3, UnstructuredGrid{3, 3, Float64, NonOriented, Nothing}}, Vector{Int64}}, AdaptedDiscreteModel{3, 3, UnstructuredDiscreteModel{3, 3, Float64, NonOriented}, UnstructuredDiscreteModel{3, 3, Float64, NonOriented}, AdaptivityGlue{Gridap.Adaptivity.RefinementGlue, 3, Vector{Vector{Int64}}, Vector{Int64}, CompressedArray{RefinementRule{Gridap.ReferenceFEs.ExtrusionPolytope{3}}, 1, Vector{RefinementRule{Gridap.ReferenceFEs.ExtrusionPolytope{3}}}, Vector{Int64}}, Table{Int64, Vector{Int64}, Vector{Int64}}, FillArrays.Fill{Bool, 1, Tuple{Base.OneTo{Int64}}}}}}, 1}, OctreeDistributedDiscreteModel{3, 3, MPIArray{Int64, 1}, GridapDistributed.GenericDistributedDiscreteModel{3, 3, MPIArray{AdaptedDiscreteModel{3, 3, UnstructuredDiscreteModel{3, 3, Float64, NonOriented}, UnstructuredDiscreteModel{3, 3, Float64, NonOriented}, AdaptivityGlue{Gridap.Adaptivity.RefinementGlue, 3, Vector{Vector{Int64}}, Vector{Int64}, CompressedArray{RefinementRule{Gridap.ReferenceFEs.ExtrusionPolytope{3}}, 1, Vector{RefinementRule{Gridap.ReferenceFEs.ExtrusionPolytope{3}}}, Vector{Int64}}, Table{Int64, Vector{Int64}, Vector{Int64}}, FillArrays.Fill{Bool, 1, Tuple{Base.OneTo{Int64}}}}}, 1}, Vector{PRange}, Nothing}, MPIArray{GridapP4est.NonConformingGlue{3, Vector{Int64}, Vector{Int64}, Vector{Vector{Tuple}}, Vector{Vector{Int64}}, Vector{Vector{Int64}}, Vector{Vector{Int64}}, Vector{Dict{Int64, Tuple{Int64, Int64, Int64}}}}, 1}, CartesianDiscreteModel{3, Float64, typeof(identity)}, Ptr{p8est_connectivity}, Ptr{p8est}}, Nothing}, cell_gids::PRange{MPIArray{LocalIndices, 1}}, cell_reffe::Gridap.ReferenceFEs.GenericLagrangianRefFE{Gridap.ReferenceFEs.GradConformity, 3}, spaces_wo_constraints::MPIArray{Gridap.FESpaces.UnconstrainedFESpace{Vector{Float64}, Gridap.FESpaces.CompressedCellConformity{Vector{Int8}}}, 1}; split_own_and_ghost::Bool, constraint::Nothing, kwargs::@Kwargs{conformity::Symbol})
# │        @ GridapP4est ~/Codes/GridapP4est.jl/src/FESpaces.jl:1061
# │     [20] _create_distributed_single_field_fe_space_with_trian_octree_model(pXest_refinement_rule_type::GridapP4est.PXestUniformRefinementRuleType, models::MPIArray{UnstructuredDiscreteModel{3, 3, Float64, NonOriented}, 1}, non_conforming_glue::MPIArray{GridapP4est.NonConformingGlue{3, Vector{Int64}, Vector{Int64}, Vector{Vector{Tuple{Int64, Int64, Int64}}}, Vector{Vector{Int64}}, Vector{Vector{Int64}}, Vector{Vector{Int64}}, Vector{Dict{Int64, Tuple{Int64, Int64, Int64}}}}, 1}, trian::DistributedTriangulation{3, 3, MPIArray{AdaptedTriangulation{3, 3, BodyFittedTriangulation{3, 3, UnstructuredDiscreteModel{3, 3, Float64, NonOriented}, GridPortion{3, 3, UnstructuredGrid{3, 3, Float64, NonOriented, Nothing}}, Vector{Int64}}, AdaptedDiscreteModel{3, 3, UnstructuredDiscreteModel{3, 3, Float64, NonOriented}, UnstructuredDiscreteModel{3, 3, Float64, NonOriented}, AdaptivityGlue{Gridap.Adaptivity.RefinementGlue, 3, Vector{Vector{Int64}}, Vector{Int64}, CompressedArray{RefinementRule{Gridap.ReferenceFEs.ExtrusionPolytope{3}}, 1, Vector{RefinementRule{Gridap.ReferenceFEs.ExtrusionPolytope{3}}}, Vector{Int64}}, Table{Int64, Vector{Int64}, Vector{Int64}}, FillArrays.Fill{Bool, 1, Tuple{Base.OneTo{Int64}}}}}}, 1}, OctreeDistributedDiscreteModel{3, 3, MPIArray{Int64, 1}, GridapDistributed.GenericDistributedDiscreteModel{3, 3, MPIArray{AdaptedDiscreteModel{3, 3, UnstructuredDiscreteModel{3, 3, Float64, NonOriented}, UnstructuredDiscreteModel{3, 3, Float64, NonOriented}, AdaptivityGlue{Gridap.Adaptivity.RefinementGlue, 3, Vector{Vector{Int64}}, Vector{Int64}, CompressedArray{RefinementRule{Gridap.ReferenceFEs.ExtrusionPolytope{3}}, 1, Vector{RefinementRule{Gridap.ReferenceFEs.ExtrusionPolytope{3}}}, Vector{Int64}}, Table{Int64, Vector{Int64}, Vector{Int64}}, FillArrays.Fill{Bool, 1, Tuple{Base.OneTo{Int64}}}}}, 1}, Vector{PRange}, Nothing}, MPIArray{GridapP4est.NonConformingGlue{3, Vector{Int64}, Vector{Int64}, Vector{Vector{Tuple}}, Vector{Vector{Int64}}, Vector{Vector{Int64}}, Vector{Vector{Int64}}, Vector{Dict{Int64, Tuple{Int64, Int64, Int64}}}}, 1}, CartesianDiscreteModel{3, Float64, typeof(identity)}, Ptr{p8est_connectivity}, Ptr{p8est}}, Nothing}, cell_gids::PRange{MPIArray{LocalIndices, 1}}, cell_reffe::MPIArray{CompressedArray{Gridap.ReferenceFEs.GenericLagrangianRefFE{Gridap.ReferenceFEs.GradConformity, 3}, 1, Vector{Gridap.ReferenceFEs.GenericLagrangianRefFE{Gridap.ReferenceFEs.GradConformity, 3}}, Vector{Int8}}, 1}; split_own_and_ghost::Bool, constraint::Nothing, kwargs::@Kwargs{conformity::Symbol})
# │        @ GridapP4est ~/Codes/GridapP4est.jl/src/FESpaces.jl:1117
# │     [21] _create_distributed_single_field_fe_space_with_trian_octree_model(_dtrian::DistributedTriangulation{3, 3, MPIArray{BodyFittedTriangulation{3, 3, UnstructuredDiscreteModel{3, 3, Float64, NonOriented}, Gridap.Geometry.GridView{3, 3, Gridap.Geometry.GridView{3, 3, UnstructuredGrid{3, 3, Float64, NonOriented, Nothing}, Vector{Int64}}, Vector{Int64}}, LazyArray{FillArrays.Fill{Reindex{IdentityVector{Int64}}, 1, Tuple{Base.OneTo{Int64}}}, Int64, 1, Tuple{LazyArray{FillArrays.Fill{Reindex{Vector{Int64}}, 1, Tuple{Base.OneTo{Int64}}}, Int64, 1, Tuple{Vector{Int64}}}}}}, 1}, OctreeDistributedDiscreteModel{3, 3, MPIArray{Int64, 1}, GridapDistributed.GenericDistributedDiscreteModel{3, 3, MPIArray{AdaptedDiscreteModel{3, 3, UnstructuredDiscreteModel{3, 3, Float64, NonOriented}, UnstructuredDiscreteModel{3, 3, Float64, NonOriented}, AdaptivityGlue{Gridap.Adaptivity.RefinementGlue, 3, Vector{Vector{Int64}}, Vector{Int64}, CompressedArray{RefinementRule{Gridap.ReferenceFEs.ExtrusionPolytope{3}}, 1, Vector{RefinementRule{Gridap.ReferenceFEs.ExtrusionPolytope{3}}}, Vector{Int64}}, Table{Int64, Vector{Int64}, Vector{Int64}}, FillArrays.Fill{Bool, 1, Tuple{Base.OneTo{Int64}}}}}, 1}, Vector{PRange}, Nothing}, MPIArray{GridapP4est.NonConformingGlue{3, Vector{Int64}, Vector{Int64}, Vector{Vector{Tuple}}, Vector{Vector{Int64}}, Vector{Vector{Int64}}, Vector{Vector{Int64}}, Vector{Dict{Int64, Tuple{Int64, Int64, Int64}}}}, 1}, CartesianDiscreteModel{3, Float64, typeof(identity)}, Ptr{p8est_connectivity}, Ptr{p8est}}, Nothing}, reffe::Tuple{Lagrangian, Tuple{DataType, Int64}, @Kwargs{}}; kwargs::@Kwargs{conformity::Symbol})
# │        @ GridapP4est ~/Codes/GridapP4est.jl/src/FESpaces.jl:1096
# │     [22] _create_distributed_single_field_fe_space_with_trian_octree_model
# │        @ ~/Codes/GridapP4est.jl/src/FESpaces.jl:1083 [inlined]
# │     [23] #FESpace#228
# │        @ ~/Codes/GridapP4est.jl/src/FESpaces.jl:1079 [inlined]
# │     [24] TestFESpace(::DistributedTriangulation{3, 3, MPIArray{BodyFittedTriangulation{3, 3, UnstructuredDiscreteModel{3, 3, Float64, NonOriented}, Gridap.Geometry.GridView{3, 3, Gridap.Geometry.GridView{3, 3, UnstructuredGrid{3, 3, Float64, NonOriented, Nothing}, Vector{Int64}}, Vector{Int64}}, LazyArray{FillArrays.Fill{Reindex{IdentityVector{Int64}}, 1, Tuple{Base.OneTo{Int64}}}, Int64, 1, Tuple{LazyArray{FillArrays.Fill{Reindex{Vector{Int64}}, 1, Tuple{Base.OneTo{Int64}}}, Int64, 1, Tuple{Vector{Int64}}}}}}, 1}, OctreeDistributedDiscreteModel{3, 3, MPIArray{Int64, 1}, GridapDistributed.GenericDistributedDiscreteModel{3, 3, MPIArray{AdaptedDiscreteModel{3, 3, UnstructuredDiscreteModel{3, 3, Float64, NonOriented}, UnstructuredDiscreteModel{3, 3, Float64, NonOriented}, AdaptivityGlue{Gridap.Adaptivity.RefinementGlue, 3, Vector{Vector{Int64}}, Vector{Int64}, CompressedArray{RefinementRule{Gridap.ReferenceFEs.ExtrusionPolytope{3}}, 1, Vector{RefinementRule{Gridap.ReferenceFEs.ExtrusionPolytope{3}}}, Vector{Int64}}, Table{Int64, Vector{Int64}, Vector{Int64}}, FillArrays.Fill{Bool, 1, Tuple{Base.OneTo{Int64}}}}}, 1}, Vector{PRange}, Nothing}, MPIArray{GridapP4est.NonConformingGlue{3, Vector{Int64}, Vector{Int64}, Vector{Vector{Tuple}}, Vector{Vector{Int64}}, Vector{Vector{Int64}}, Vector{Vector{Int64}}, Vector{Dict{Int64, Tuple{Int64, Int64, Int64}}}}, 1}, CartesianDiscreteModel{3, Float64, typeof(identity)}, Ptr{p8est_connectivity}, Ptr{p8est}}, Nothing}, ::Vararg{Any}; kwargs::@Kwargs{conformity::Symbol})
# │        @ Gridap.FESpaces ~/.julia/packages/Gridap/ZgUBG/src/FESpaces/FESpaceFactories.jl:5
# │     [25] (::var"#mwe_sphere_3d_hanging_glue_on_portion##0#mwe_sphere_3d_hanging_glue_on_portion##1"{Int64, Int64, MPIArray{Int64, 1}})()
# │        @ Main ~/Codes/GridapP4est.jl/test/issue_108.jl:76
# │     [26] #with#2
# │        @ ~/Codes/GridapP4est.jl/src/Environment.jl:41 [inlined]
# │     [27] with
# │        @ ~/Codes/GridapP4est.jl/src/Environment.jl:39 [inlined]
# │     [28] mwe_sphere_3d_hanging_glue_on_portion(distribute::PartitionedArrays.var"#60#61"{MPI.Comm, Bool}, nparts::Int64, n::Int64)
# │        @ Main ~/Codes/GridapP4est.jl/test/issue_108.jl:25
# │     [29] #6
# │        @ ~/Codes/GridapP4est.jl/test/issue_108.jl:94 [inlined]
# │     [30] with_mpi(f::var"#6#7"; comm::MPI.Comm, duplicate_comm::Bool)
# │        @ PartitionedArrays ~/.julia/packages/PartitionedArrays/MVmxR/src/mpi_array.jl:73
# │     [31] with_mpi(f::Function)
# │        @ PartitionedArrays ~/.julia/packages/PartitionedArrays/MVmxR/src/mpi_array.jl:64
# │     [32] top-level scope
# │        @ ~/Codes/GridapP4est.jl/test/issue_108.jl:93
# │     [33] include(mod::Module, _path::String)
# │        @ Base ./Base.jl:306
# │     [34] exec_options(opts::Base.JLOptions)
# │        @ Base ./client.jl:317
# │     [35] _start()
# │        @ Base ./client.jl:550
# └ @ PartitionedArrays ~/.julia/packages/PartitionedArrays/MVmxR/src/mpi_array.jl:75
# Abort(1) on node 0 (rank 0 in comm 0): application called MPI_Abort(MPI_COMM_WORLD, 1) - process 0
   
    e = u - uh
    el2 = sqrt(∑( ∫( e⋅e )dΩ1 ))
    println("L2 error: $el2")

    # writevtk(Ω1, joinpath(outdir, "Ω1"), cellfields=["uh"=>uh,"u_ex"=>u,"error"=>e])

    return nothing

    end # GridapP4est.with
end


with_mpi() do distribute
    mwe_sphere_3d_hanging_glue_on_portion(distribute,3,3)
    println("finished")
end