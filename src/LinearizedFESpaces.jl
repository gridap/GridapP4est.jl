
function _get_order(reffe::Tuple{<:Lagrangian,Any,Any})
    reffe[2][2]
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
        ref_model,glue=Gridap.adapt(ref_model,ref_coarse_flags)
        #writevtk(ref_model,"ref_model")
    end
    Gridap.FESpace(ref_model,lin_reffe; kwargs...), ref_model
end