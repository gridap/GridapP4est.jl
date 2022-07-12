struct FESpaceHierarchyLevel{A,B}
  level        :: Int
  fe_space     :: A
  fe_space_red :: B
end

function get_level_fe_space(a::FESpaceHierarchyLevel{A,Nothing}) where {A}
  a.fe_space
end

function get_level_fe_space(a::FESpaceHierarchyLevel{A,B}) where {A,B}
  a.fe_space_red
end

function get_level_fe_space_before_redist(a::FESpaceHierarchyLevel)
  a.fe_space
end

function Gridap.FESpaces.TestFESpace(
      mh::ModelHierarchyLevel{A,B,Nothing},args...;kwargs...) where {A,B}
  Vh = TestFESpace(mh.model.dmodel,args...;kwargs...)
  FESpaceHierarchyLevel(mh.level,Vh,nothing)
end

function Gridap.FESpaces.TestFESpace(
      mh::ModelHierarchyLevel{A,B,C},args...;kwargs...) where {A,B,C}
  Vh     = TestFESpace(mh.model.dmodel,args...;kwargs...)
  Vh_red = TestFESpace(mh.model_red.dmodel,args...;kwargs...)
  FESpaceHierarchyLevel(mh.level,Vh,Vh_red)
end

function Gridap.FESpaces.TrialFESpace(u,a::FESpaceHierarchyLevel{A,Nothing}) where {A}
  Uh = TrialFESpace(u,a.fe_space)
  FESpaceHierarchyLevel(a.level,Uh,nothing)
end

function Gridap.FESpaces.TrialFESpace(u,a::FESpaceHierarchyLevel{A,B}) where {A,B}
  Uh     = TrialFESpace(u,a.fe_space)
  Uh_red = TrialFESpace(u,a.fe_space_red)
  FESpaceHierarchyLevel(a.level,Uh,Uh_red)
end
