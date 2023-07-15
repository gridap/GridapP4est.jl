function Gridap.CellData.change_domain(a::CellField,target_trian::Triangulation,target_domain::DomainStyle)
  strian=get_triangulation(a) 
  if (strian===target_trian)
    change_domain(a,DomainStyle(a),target_domain)
  else
    change_domain(a,get_triangulation(a),DomainStyle(a),target_trian,target_domain)
  end
end