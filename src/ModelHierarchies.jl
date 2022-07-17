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
get_level(a::ModelHierarchy,level::Integer) = a.levels[level]

get_level_model(a::ModelHierarchy,level::Integer)=get_level_model(get_level(a,level))
get_level_model(a::ModelHierarchyLevel{A,B,Nothing}) where {A,B}=a.model
get_level_model(a::ModelHierarchyLevel{A,B,C}) where {A,B,C}=a.model_red

get_level_model_before_redist(a::ModelHierarchy,level::Integer)=
     get_level_model_before_redist(get_level(a,level))
get_level_model_before_redist(a::ModelHierarchyLevel) where {A,B}=a.model

# Implement support for num_refs_x_level? (future work)
function ModelHierarchy(parts,cmodel,num_procs_x_level; num_refs_x_level=nothing)
  num_levels  = length(num_procs_x_level)
  level_parts = generate_level_parts(parts,num_procs_x_level)

  model=OctreeDistributedDiscreteModel(level_parts[num_levels],cmodel)
  meshes=Vector{ModelHierarchyLevel}(undef,num_levels)

  meshes[num_levels]=ModelHierarchyLevel(num_levels,model,nothing,nothing,nothing)

  for i=num_levels-1:-1:1
    modelH=get_level_model(meshes[i+1])
    if (num_procs_x_level[i]!=num_procs_x_level[i+1])
      # meshes[i+1].model is distributed among P processors
      # model_ref is distributed among Q processors, with P!=Q
      model_ref,ref_glue=refine(modelH,level_parts[i])
      model_red,red_glue=redistribute(model_ref)
    else
      model_ref,ref_glue=refine(modelH)
      model_red,red_glue=nothing,nothing
    end
    meshes[i]=ModelHierarchyLevel(i,model_ref,ref_glue,model_red,red_glue)
  end
  ModelHierarchy(level_parts,meshes)
end
