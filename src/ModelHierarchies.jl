# 1. redistributed and red_glue might be of type Nothing
#    whenever there is no redistribution in a given level
# 2. ref_glue is of type Nothing for the coarsest model
struct ModelHierarchyLevel{A,B,C,D}
  level::Int
  model::A
  ref_glue::B
  model_red::C
  red_glue::D
end

function model_hierarchy_level_free!(a::ModelHierarchyLevel{A,B,Nothing,Nothing}) where {A,B}
  octree_distributed_discrete_model_free!(a.model)
end
function model_hierarchy_level_free!(a::ModelHierarchyLevel{A,B,C,D}) where {A,B,C,D}
  octree_distributed_discrete_model_free!(a.model)
  octree_distributed_discrete_model_free!(a.model_red)
end

struct ModelHierarchy
  level_parts::Vector{PArrays.MPIData}
  levels::Vector{ModelHierarchyLevel}
end

function model_hierarchy_free!(a::ModelHierarchy)
  for level in a.levels
    model_hierarchy_level_free!(level)
  end
end


num_levels(a::ModelHierarchy)= length(a.levels)
get_level(a::ModelHierarchy,level) = a.levels[level]

# Implement support for num_refs_x_level? (future work)
function ModelHierarchy(parts,cmodel,num_procs_x_level; num_refs_x_level=nothing)
  num_levels  = length(num_procs_x_level)
  level_parts = generate_level_parts(parts,num_procs_x_level)

  model=OctreeDistributedDiscreteModel(level_parts[num_levels],cmodel)
  meshes=Vector{ModelHierarchyLevel}(undef,num_levels)

  meshes[num_levels]=ModelHierarchyLevel(num_levels,model,nothing,nothing,nothing)

  for i=num_levels-1:-1:1
    model_ref,ref_glue=refine(meshes[i+1].model)
    if (num_procs_x_level[i]!=num_procs_x_level[i+1])
      model_red,red_glue=redistribute(model_ref,level_parts[i])
    else
      model_red,red_glue=nothing,nothing
    end
    meshes[i]=ModelHierarchyLevel(i,model_ref,ref_glue,model_red,red_glue)
  end
  ModelHierarchy(level_parts,meshes)
end
