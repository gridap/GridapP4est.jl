function _create_ref_rule(Dc,ranks::MPIArray, num_uniform_refs)
    if Dc==2
      coarse_model=CartesianDiscreteModel((0,1,0,1),(1,1))
    else
      @assert Dc==3
      coarse_model=CartesianDiscreteModel((0,1,0,1,0,1),(1,1,1))
    end   
    new_comm=MPI.Comm_split(ranks.comm,MPI.Comm_rank(ranks.comm),0)
    new_ranks=MPIArray(1,new_comm,(1,))
    model_ref=OctreeDistributedDiscreteModel(new_ranks,coarse_model,num_uniform_refs)
    ref_rules=map(local_views(model_ref.dmodel)) do model
      Gridap.Adaptivity.RefinementRule(Gridap.Adaptivity.GenericRefinement(),
                                       Dc==2 ? QUAD : HEX, 
                                       get_grid(model))
    end
    ref_rules.item_ref[]
end

function _get_order(reffe::Tuple{<:Lagrangian,Any,Any})
    reffe[2][2]
end

function _create_adaptivity_glue(model::OctreeDistributedDiscreteModel{Dc},
                                 ref_model::OctreeDistributedDiscreteModel{Dc},
                                 num_uniform_refinements) where {Dc}
    ref_rule=_create_ref_rule(Dc,model.parts,num_uniform_refinements)
    num_children=Gridap.Adaptivity.num_subcells(ref_rule)
    cell_gids_model     = get_cell_gids(model.dmodel)
    cell_gids_ref_model = get_cell_gids(ref_model.dmodel)

    n2o_cell_map,n2o_cell_to_child_id=map(partition(cell_gids_model),partition(cell_gids_ref_model)) do model_partition, ref_model_partition
        num_local_cells_model=local_length(model_partition)
        num_owned_cells_model=own_length(model_partition)
        num_local_cells_ref_model=local_length(ref_model_partition)
        n2o_cell_map=Vector{Int}(undef,num_local_cells_ref_model)
        n2o_cell_to_child_id=Vector{Int}(undef,num_local_cells_ref_model)
        current=1
        for i=1:num_owned_cells_model
            for j=1:num_children
               n2o_cell_map[current]=i
               n2o_cell_to_child_id[current]=j
               current+=1
            end
        end
        n2o_cell_map,n2o_cell_to_child_id
    end |> tuple_of_arrays

    cache = fetch_vector_ghost_values_cache(n2o_cell_map,partition(cell_gids_ref_model))
    fetch_vector_ghost_values!(n2o_cell_map,cache) |> wait
    fetch_vector_ghost_values!(n2o_cell_to_child_id,cache) |> wait

    adaptivity_glue=map(n2o_cell_map,n2o_cell_to_child_id) do n2o_cell_map, n2o_cell_to_child_id
        n2o_faces_map = [(d==Dc) ? n2o_cell_map : Int[] for d in 0:Dc]
        AdaptivityGlue(n2o_faces_map,n2o_cell_to_child_id,ref_rule)
    end
end 

function _setup_one_level_refined_octree_model(model::OctreeDistributedDiscreteModel{Dc,Dp},
                                               fmodel::OctreeDistributedDiscreteModel{Dc,Dp},
                                               adaptivity_glue) where {Dc,Dp}
    @assert model.parts === fmodel.parts
    adaptive_models = map(local_views(model),
                          local_views(fmodel),
                          adaptivity_glue) do model, fmodel, glue 
      Gridap.Adaptivity.AdaptedDiscreteModel(fmodel.model,model,glue)
    end
    new_fmodel = GridapDistributed.GenericDistributedDiscreteModel(adaptive_models,get_cell_gids(fmodel))
    OctreeDistributedDiscreteModel(Dc,Dp,
                        model.parts,
                        new_fmodel,
                        fmodel.non_conforming_glue,
                        model.coarse_model,
                        model.ptr_pXest_connectivity,
                        pXest_copy(fmodel.pXest_type,fmodel.ptr_pXest),
                        fmodel.pXest_type,
                        fmodel.pXest_refinement_rule_type,
                        false,
                        fmodel)
end 

function Gridap.LinearizedFESpace(model::OctreeDistributedDiscreteModel{Dc}, 
                                  reffe::Tuple{Gridap.ReferenceFEs.Lagrangian,Any,Any}; 
                                  kwargs...) where {Dc}
    lin_reffe=Gridap._linearize_reffe(reffe)
    order=_get_order(reffe)
    @assert floor(log2(order)) == ceil(log2(order)) "The order of the Lagrangian reference FE must be a power of 2"
    num_refinements=Int(log2(order))
    ref_model=model
    for i=1:num_refinements
        cell_gids=get_cell_gids(ref_model.dmodel)
        ref_coarse_flags=map(partition(cell_gids)) do indices
            flags=Vector{Cint}(undef,local_length(indices))
            flags.=refine_flag
        end
        ref_model,glue=Gridap.Adaptivity.adapt(ref_model,ref_coarse_flags)
    end
    adaptivity_glue=_create_adaptivity_glue(model,ref_model,num_refinements)
    one_level_ref_model=_setup_one_level_refined_octree_model(model,ref_model,adaptivity_glue)
    Gridap.FESpace(one_level_ref_model,lin_reffe; kwargs...), one_level_ref_model
end