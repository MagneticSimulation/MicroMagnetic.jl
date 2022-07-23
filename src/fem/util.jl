function kronecker_delta(x::Int64, y::Int64)
  if x == y 
    return 1
  end
  return 0
end

#
# helper function, compute the det of the given three vectors
# for equation (9), equals 6*a1*Ve
#
function det3(x::Array{Float64,1}, y::Array{Float64,1}, z::Array{Float64,1})
  d = x[1] * y[2] * z[3] + x[2] * y[3] * z[1] + x[3] * y[1] * z[2] -
      x[1] * y[3] * z[2] - x[2] * y[1] * z[3] - x[3] * y[2] * z[1];
  return d;
end

function norm2(X::Array{Float64,1})
  v = 0.0
  for x in X
      v += x*x
  end
  return sqrt(v)
end

function dot3(x::Array{Float64,1}, y::Array{Float64,1})
  return x[1]*y[1] + x[2]*y[2] + x[3]*y[3] 
end

function triangle_area(x1::Array{Float64,1}, x2::Array{Float64,1}, x3::Array{Float64,1})
  a = x2 .- x1
  b = x3 .- x2
  c = x1 .- x3
  s1 = norm2(a);
  s2 = norm2(b);
  s3 = norm2(c);
  s = (s1 + s2 + s3) / 2.0;
  area = sqrt(s * (s - s1) * (s - s2) * (s - s3))
  return area
end


function cross3(x::Array{Float64,1}, y::Array{Float64,1})
  return [-x[3]*y[2] + x[2]*y[3], x[3]*y[1] - x[1]*y[3], -x[2]*y[1] + x[1]*y[2]]
end

#compute the volume of a tetrahedron
function tetrahedron_volume(x1::Array{Float64,1}, x2::Array{Float64,1}, x3::Array{Float64,1}, x4::Array{Float64,1})
  v = det3(x2, x3, x4) - det3(x1, x3, x4) + det3(x1, x2, x4) - det3(x1, x2, x3);
  return v / 6.0;
end


#
#equation (10), (11) and (12), note the factor of 6*Ve
#compute the coefficients b c and d of a tetrahedron
#return (b1,c1,d1)
#
function bcd_negative(x2::Array{Float64,1}, x3::Array{Float64,1}, x4::Array{Float64,1})
    
    b1 = x2[3] * x3[2] - x2[2] * x3[3] - x2[3] * x4[2] + x3[3] * x4[2] + x2[2] * x4[3] - x3[2] * x4[3];

    c1 = -x2[3] * x3[1] + x2[1] * x3[3] + x2[3] * x4[1] - x3[3] * x4[1] - x2[1] * x4[3] + x3[1] * x4[3];

    d1 = x2[2] * x3[1] - x2[1] * x3[2] - x2[2] * x4[1] + x3[2] * x4[1] + x2[1] * x4[2] - x3[1] * x4[2];

    return (b1, c1, d1)
end


#
# compute alpha \times k \times beta, where alpha,beta=1,2,3,4; k=1,2,3;
# return b1 c1 d1 b2 c2 d2 b3 c3 d3 b4 c4 d4
# 
function compute_bcd12(x1::Array{Float64,1}, x2::Array{Float64,1}, x3::Array{Float64,1}, x4::Array{Float64,1})

  J = 1.0/(6*tetrahedron_volume(x1, x2, x3, x4))

  (b1, c1, d1) = bcd_negative(x2, x3, x4) .* J
  (b2, c2, d2) = bcd_negative(x1, x3, x4) .* (-J)
  (b3, c3, d3) = bcd_negative(x2, x1, x4) .* (-J)
  (b4, c4, d4) = bcd_negative(x2, x3, x1) .* (-J)

  return (b1,c1,d1,b2,c2,d2,b3,c3,d3,b4,c4,d4)
end

# compute the mesh volumes of given mesh
function compute_mesh_volume!(mesh::FEMesh)

  volumes = zeros(mesh.number_cells)

  for i = 1:mesh.number_cells

    k1 = mesh.cell_verts[1, i]
    k2 = mesh.cell_verts[2, i]
    k3 = mesh.cell_verts[3, i]
    k4 = mesh.cell_verts[4, i]

    v = tetrahedron_volume(mesh.coordinates[:, k1], 
                           mesh.coordinates[:, k2],
                           mesh.coordinates[:, k3],
                           mesh.coordinates[:, k4])
    
    volumes[i]  =  abs(v)

  end

  mesh.volumes = volumes

end

function compute_L_inv_neg!(mesh::FEMesh)
  nodal_volumes = zeros(mesh.number_nodes)
  cv = mesh.cell_verts
  for i = 1:mesh.number_cells, k = 1:4
    nodal_volumes[cv[k,i]] += mesh.volumes[i]/4.0  
  end

  L_inv_neg = zeros(3*mesh.number_nodes)

  for i = 0:mesh.number_nodes-1
    L_inv_neg[3*i+1] = -1.0/nodal_volumes[i+1]
    L_inv_neg[3*i+2] = -1.0/nodal_volumes[i+1]
    L_inv_neg[3*i+3] = -1.0/nodal_volumes[i+1]
  end

  mesh.L_inv_neg = L_inv_neg

end

function compute_L_Ms!(L_mu::Array{Float64,1}, mesh::FEMesh, Ms::Array{Float64,1})
  nodal_Ms = zeros(mesh.number_nodes)
  cv = mesh.cell_verts

  for k = 1:4, i = 1:mesh.number_cells
    nodal_Ms[cv[k,i]] += mesh.volumes[i]*Ms[i]/4.0  
  end

  for i = 0:mesh.number_nodes-1
    L_mu[3*i+1] = nodal_Ms[i+1]
    L_mu[3*i+2] = nodal_Ms[i+1]
    L_mu[3*i+3] = nodal_Ms[i+1]
  end

end

#areas is an array
#property is an array with the same length of areas and contain
#the property information
#target is used to filter areas according to property
function max_id_with_filter(areas, property, target)
  max_value, id = minimum(areas), -1
  for i = 1:length(areas)
    if areas[i] >= max_value && property[i] == target
      max_value = areas[i]
      id = i
    end
  end
  return id
end

function compute_normal(x1::Array{Float64,1}, x2::Array{Float64,1}, x3::Array{Float64,1})
  v1 = x2 .- x1
  v2 = x3 .- x1
  v = cross3(v1, v2)
  return v/norm2(v)
end

function extract_normal_axes_by_maximum_area(mesh::FEMesh)

  N_surfaces = length(Set(mesh.surface_ids)) #the total surfaces number of the defined geometry
  N_materials = length(Set(mesh.material_ids)) #the total materials number

  surface_normals = zeros(3, N_surfaces)
  surface_node = zeros(Int64, N_surfaces) # we only need one node for each surface

  areas = zeros(N_surfaces)
  #compute the total areas of each surfaces
  for f = 1:mesh.number_faces_bnd
    i = mesh.face_verts[1, f]
    j = mesh.face_verts[2, f]
    k = mesh.face_verts[3, f]

    v1 = mesh.coordinates[:, i]
    v2 = mesh.coordinates[:, j]
    v3 = mesh.coordinates[:, k]

    id = mesh.surface_ids[f]
    if surface_node[id] == 0  # we calculate the surface normal and store one node to surface_node array
      surface_node[id] = i
      surface_normals[:, id] .= compute_normal(v1, v2, v3)
    end 
    areas[id] += triangle_area(v1, v2, v3)
  end

    #mapping for materials and nodes
    materials_dict = Dict{Int64, Set}()
    for c = 1:mesh.number_cells
      mid = mesh.material_ids[c]
      if !haskey(materials_dict, mid)
        materials_dict[mid] = Set{Int64}()
      end
      for i in mesh.cell_verts[:, c]
          push!(materials_dict[mid], i)
      end
    end

  #we need to group the areas according to material ids.
  #we assume that each individual geometry does not connect to other geometries
  #so each surface should only belongs to one material
  map_surface2material = zeros(Int64, N_surfaces)
  for i = 1:N_surfaces
    for j = 1:length(materials_dict)
      if surface_node[i] in materials_dict[j]
        map_surface2material[i] = j
        break
      end
    end
  end


  #finally, we find the maximum area for each material and return the corresponding normals.
  max_normals = zeros(3, N_materials)
  for i = 1:N_materials
    id = max_id_with_filter(areas, map_surface2material, i)
    max_normals[:,i] .= surface_normals[:, id]
  end

  return max_normals

end

