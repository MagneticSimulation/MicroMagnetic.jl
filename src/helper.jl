function  init_scalar!(v::Array{T, 1}, mesh::Mesh, init::Number) where {T<:AbstractFloat}
    for i = 1:mesh.nxyz
        v[i] = init
    end
    return nothing
end

function  init_scalar!(v::Array{T, 1}, mesh::Mesh, init::Array{T, 1}) where {T<:AbstractFloat}
    v .= init
    return nothing
end

function init_scalar!(v::Array{T, 1}, mesh::Mesh, init_fun::Function) where {T<:AbstractFloat}
    mesh = mesh
    for k = 1:mesh.nz, j = 1:mesh.ny, i = 1:mesh.nx
        id = index(i, j, k, mesh.nx, mesh.ny, mesh.nz)
        v[id] = init_fun(i, j, k, mesh.dx, mesh.dy, mesh.dz)
    end
    return true
end

function init_vector!(v::Array{T, 1}, mesh::Mesh, init::Function) where {T<:AbstractFloat}
  nxyz = mesh.nx*mesh.ny*mesh.nz
  dx,dy,dz = mesh.dx, mesh.dy, mesh.dz
  b = reshape(v, 3, nxyz)
  for i = 1:mesh.nx
    for j = 1:mesh.ny
      for k = 1:mesh.nz
        id = index(i, j, k, mesh.nx, mesh.ny, mesh.nz)
        b[:, id] .=  init(i,j,k,dx,dy,dz)
      end
    end
  end
  if NaN in v
      error("NaN is given by the input function.")
  end
  return nothing
end

function init_vector!(v::Array{T, 1}, mesh::CubicMesh, init::Function) where {T<:AbstractFloat}
  nxyz = mesh.nx*mesh.ny*mesh.nz
  b = reshape(v, 3, nxyz)
  for i = 1:mesh.nx
    for j = 1:mesh.ny
      for k = 1:mesh.nz
        id = index(i, j, k, mesh.nx, mesh.ny, mesh.nz)
        b[:, id] .=  init(i,j,k)
      end
    end
  end
end

function init_vector!(v::Array{T, 1}, mesh::Mesh, init::Tuple{Real,Real,Real}) where {T<:AbstractFloat}
  nxyz = mesh.nx*mesh.ny*mesh.nz
  b = reshape(v, 3, nxyz)
  b[1, :] .= init[1]
  b[2, :] .= init[2]
  b[3, :] .= init[3]
  return nothing
end

function init_vector!(v::Array{T, 1}, mesh::Mesh, init::Array{T, 1}) where {T<:AbstractFloat}
  v .= init
  return nothing
end

function normalise(a::Array{T, 1}, N::Int64) where {T<:AbstractFloat}
  for i = 0:N-1
      j = 3*i+1
      length = sqrt(a[j]*a[j]+a[j+1]*a[j+1]+a[j+2]*a[j+2])
      if length > 0
        length = 1.0/length
        a[j] *= length
        a[j+1] *= length
        a[j+2] *= length
      end
    end
   return nothing
end

function error_length_m(a::Array{Float64, 1}, N::Int64)
  maxlength = 0.0
  minlength = 1.0
  for i = 0:N-1
    j = 3*i+1
    length = sqrt(a[j]*a[j]+a[j+1]*a[j+1]+a[j+2]*a[j+2])
    if length > maxlength
      maxlength = length
    end
    if length < minlength
      minlength = length
    end
  end
  return maxlength-minlength
end

function omega_to_spin(omega::Array{Float64, 1}, spin::Array{Float64, 1}, spin_next::Array{Float64, 1}, N::Int64)
  #compute Cay(Omega).m where Cay(Omega) = (I - 1/2 Omega)^-1 (I + 1/2 Omega)
  #where Omega = Skew[w1, w2, w3] = {{0, -w3, w2}, {w3, 0, -w1}, {-w2, w1, 0}}
  for i = 0:N-1
    j = 3*i + 1
    w1 = omega[j]*0.5
    w2 = omega[j+1]*0.5
    w3 = omega[j+2]*0.5
    m1 = spin[j]
    m2 = spin[j+1]
    m3 = spin[j+2]
    r = 1 + w1*w1 + w2*w2 + w3*w3
    a11 = 1 + w1*w1 - w2*w2 - w3*w3
    a12 = 2*(w1*w2 - w3)
    a13 = 2*(w2 + w1*w3)
    a21 = 2*(w1*w2 + w3)
    a22 = 1 - w1*w1 + w2*w2 - w3*w3
    a23 = -2*(w1-w2*w3)
    a31 = 2*(-w2+w1*w3)
    a32 = 2*(w1+w2*w3)
    a33 = 1 - w1*w1 - w2*w2 + w3*w3
    spin_next[j] = (a11*m1 + a12*m2 + a13*m3)/r
    spin_next[j+1] = (a21*m1 + a22*m2 + a23*m3)/r
    spin_next[j+2] = (a31*m1 + a32*m2 + a33*m3)/r
  end
end
function compute_error2(error::Array{Float64,1}, N::Int64)
  norm = 0.0
  for i=1:N
    norm += error[i]^2
  end
  return sqrt(norm/N)
end

function compute_dmdt(m1::Array{Float64, 1}, m2::Array{Float64, 1}, N::Int64, dt::Float64)
  max_dmdt = 0.0
  for i = 0:N-1
    j = 3*i + 1
    mx = m1[j] - m2[j]
    my = m1[j+1] - m2[j+1]
    mz = m1[j+2] - m2[j+2]
    dmdt = sqrt(mx*mx + my*my + mz*mz)
    if dmdt > max_dmdt
      max_dmdt = dmdt
    end
  end
  return max_dmdt/dt
end

function compute_dm_step(m1::Array{Float64, 1}, m2::Array{Float64, 1}, N::Int64)
  max_dm = 0.0
  for i = 0:N-1
    j = 3*i + 1
    mx = m1[j] - m2[j]
    my = m1[j+1] - m2[j+1]
    mz = m1[j+2] - m2[j+2]
    dmdt = sqrt(mx*mx + my*my + mz*mz)
    if dmdt > max_dm
      max_dm = dmdt
    end
  end
  return max_dm
end


function compute_skyrmion_number(v::Array{T, 1}, m::Array{T, 1}, mesh::Mesh) where {T<:AbstractFloat}
    nx,ny,nz = mesh.nx, mesh.ny, mesh.nz
    for k = 1:nz, j = 1:ny, i=1:nx
        id = index(i, j, k, nx, ny, nz)
        mx,my,mz = m[3*id-2],m[3*id-1],m[3*id]
        sx1,sy1,sz1 = T(0),T(0),T(0)
        sx2,sy2,sz2 = T(0),T(0),T(0)
        id1 = 3*_x_minus_one(i, id, nx, ny, nz, mesh.xperiodic)
        id2 = 3*_y_minus_one(j, id, nx, ny, nz, mesh.yperiodic)
        v[id] = 0
        if id1>0 && id2>0
            sx1,sy1,sz1 = m[id1-2],m[id1-1],m[id1]
            sx2,sy2,sz2 = m[id2-2],m[id2-1],m[id2]
            v[id] += Berg_Omega(sx2, sy2, sz2, mx, my, mz, sx1, sy1, sz1)
        end

        id1 = 3*_x_plus_one(i, id, nx, ny, nz, mesh.xperiodic)
        id2 = 3*_y_plus_one(j, id, nx, ny, nz, mesh.yperiodic)
        if id1>0 && id2>0
            sx1,sy1,sz1 = m[id1-2],m[id1-1],m[id1]
            sx2,sy2,sz2 = m[id2-2],m[id2-1],m[id2]
            v[id] += Berg_Omega(sx2, sy2, sz2, mx, my, mz, sx1, sy1, sz1)
        end
        v[id] /= (4*pi);
    end
    return nothing
end

function compute_winding_number_xy(m::Array{T, 1}, mesh::Mesh; k=1) where {T<:AbstractFloat}
    nx,ny,nz = mesh.nx, mesh.ny, mesh.nz
    v = 0
    for j = 1:ny, i=1:nx
        id = index(i, j, k, nx, ny, nz)
        mx,my,mz = m[3*id-2],m[3*id-1],m[3*id]
        sx1,sy1,sz1 = T(0),T(0),T(0)
        sx2,sy2,sz2 = T(0),T(0),T(0)
        id1 = 3*_x_minus_one(i, id, nx, ny, nz, mesh.xperiodic)
        id2 = 3*_y_minus_one(j, id, nx, ny, nz, mesh.yperiodic)
        if id1>0 && id2>0
            sx1,sy1,sz1 = m[id1-2],m[id1-1],m[id1]
            sx2,sy2,sz2 = m[id2-2],m[id2-1],m[id2]
            v += Berg_Omega(sx2, sy2, sz2, mx, my, mz, sx1, sy1, sz1)
        end

        id1 = 3*_x_plus_one(i, id, nx, ny, nz, mesh.xperiodic)
        id2 = 3*_y_plus_one(j, id, nx, ny, nz, mesh.yperiodic)
        if id1>0 && id2>0
            sx1,sy1,sz1 = m[id1-2],m[id1-1],m[id1]
            sx2,sy2,sz2 = m[id2-2],m[id2-1],m[id2]
            v += Berg_Omega(sx2, sy2, sz2, mx, my, mz, sx1, sy1, sz1)
        end

    end
    return v/(4*pi);
end


function compute_winding_number_yz(m::Array{T, 1}, mesh::Mesh; i=1) where {T<:AbstractFloat}
    nx,ny,nz = mesh.nx, mesh.ny, mesh.nz
    v = 0
    for k=1:nz, j=1:ny
        id = index(i, j, k, nx, ny, nz)
        mx,my,mz = m[3*id-2],m[3*id-1],m[3*id]
        sx1,sy1,sz1 = T(0),T(0),T(0)
        sx2,sy2,sz2 = T(0),T(0),T(0)
        id1 = 3*_y_minus_one(j, id, nx, ny, nz, mesh.yperiodic)
        id2 = 3*_z_minus_one(k, id, nx, ny, nz, mesh.zperiodic)
        if id1>0 && id2>0
            sx1,sy1,sz1 = m[id1-2],m[id1-1],m[id1]
            sx2,sy2,sz2 = m[id2-2],m[id2-1],m[id2]
            v += Berg_Omega(sx2, sy2, sz2, mx, my, mz, sx1, sy1, sz1)
        end

        id1 = 3*_y_plus_one(j, id, nx, ny, nz, mesh.yperiodic)
        id2 = 3*_z_plus_one(k, id, nx, ny, nz, mesh.zperiodic)
        if id1>0 && id2>0
            sx1,sy1,sz1 = m[id1-2],m[id1-1],m[id1]
            sx2,sy2,sz2 = m[id2-2],m[id2-1],m[id2]
            v += Berg_Omega(sx2, sy2, sz2, mx, my, mz, sx1, sy1, sz1)
        end

    end
    return v/(4*pi);
end


function compute_winding_number_xz(m::Array{T, 1}, mesh::Mesh; j=1) where {T<:AbstractFloat}
    nx,ny,nz = mesh.nx, mesh.ny, mesh.nz
    v = 0
    for k=1:nz, i=1:nx
        id = index(i, j, k, nx, ny, nz)
        mx,my,mz = m[3*id-2],m[3*id-1],m[3*id]
        sx1,sy1,sz1 = T(0),T(0),T(0)
        sx2,sy2,sz2 = T(0),T(0),T(0)
        id1 = 3*_z_minus_one(k, id, nx, ny, nz, mesh.zperiodic)
        id2 = 3*_x_minus_one(i, id, nx, ny, nz, mesh.xperiodic)
        if id1>0 && id2>0
            sx1,sy1,sz1 = m[id1-2],m[id1-1],m[id1]
            sx2,sy2,sz2 = m[id2-2],m[id2-1],m[id2]
            v += Berg_Omega(sx2, sy2, sz2, mx, my, mz, sx1, sy1, sz1)
        end

        id1 = 3*_z_plus_one(k, id, nx, ny, nz, mesh.zperiodic)
        id2 = 3*_x_plus_one(i, id, nx, ny, nz, mesh.xperiodic)
        if id1>0 && id2>0
            sx1,sy1,sz1 = m[id1-2],m[id1-1],m[id1]
            sx2,sy2,sz2 = m[id2-2],m[id2-1],m[id2]
            v += Berg_Omega(sx2, sy2, sz2, mx, my, mz, sx1, sy1, sz1)
        end

    end
    return v/(4*pi);
end


#for a cubic mesh, each cube has six faces (12 triangles)
#Here we can construct an Octahedron for each site, so we have 8 triangles, we only need two?
function compute_winding_number_3d(v::Array{T, 1}, m::Array{T, 1}, mesh::Mesh) where {T<:AbstractFloat}
    nx,ny,nz = mesh.nx, mesh.ny, mesh.nz
    for k = 1:nz, j = 1:ny, i=1:nx
        id = index(i, j, k, nx, ny, nz)
        v[id] = 0

        mx,my,mz = m[3*id-2],m[3*id-1],m[3*id]
        sx1,sy1,sz1 = T(0),T(0),T(0)
        sx2,sy2,sz2 = T(0),T(0),T(0)
        sx3,sy3,sz3 = T(0),T(0),T(0)

        #(-x, -y, -z) triangle
        id1 = 3*_x_minus_one(i, id, nx, ny, nz, mesh.xperiodic)
        id2 = 3*_y_minus_one(j, id, nx, ny, nz, mesh.yperiodic)
        id3 = 3*_z_minus_one(k, id, nx, ny, nz, mesh.zperiodic)
        if id1>0 && id2>0 && id3>0
            sx1,sy1,sz1 = m[id1-2],m[id1-1],m[id1]
            sx2,sy2,sz2 = m[id2-2],m[id2-1],m[id2]
            sx3,sy3,sz3 = m[id3-2],m[id3-1],m[id3]
            v[id] += Berg_Omega(sx1, sy1, sz1, sx2, sy2, sz2, sx3, sy3, sz3)
        end

        #(+x, +y, +z) triangle
        id1 = 3*_x_plus_one(i, id, nx, ny, nz, mesh.xperiodic)
        id2 = 3*_y_plus_one(j, id, nx, ny, nz, mesh.yperiodic)
        id3 = 3*_z_plus_one(k, id, nx, ny, nz, mesh.zperiodic)
        if id1>0 && id2>0 && id3>0
            sx1,sy1,sz1 = m[id1-2],m[id1-1],m[id1]
            sx2,sy2,sz2 = m[id2-2],m[id2-1],m[id2]
            sx3,sy3,sz3 = m[id3-2],m[id3-1],m[id3]
            v[id] -= Berg_Omega(sx1, sy1, sz1, sx2, sy2, sz2, sx3, sy3, sz3)
        end

        v[id] /= (2*pi);
    end
    return nothing
end

function compute_winding_number_3d(m::Array{T, 1}, mesh::Mesh) where {T<:AbstractFloat}
    nx,ny,nz = mesh.nx, mesh.ny, mesh.nz
    v = zeros(T, 3*ny*ny*nz)
    ompute_winding_number_3d(v, m, mesh)
    return sum(v)
end

function winding_number_3d(m::Array{T, 1}, mesh::Mesh) where {T<:AbstractFloat}
    nx,ny,nz = mesh.nx, mesh.ny, mesh.nz
    v = 0
    v -= compute_winding_number_xy(m, mesh, k=1)
    v += compute_winding_number_xy(m, mesh, k=nz)

    v -= compute_winding_number_yz(m, mesh, i=1)
    v += compute_winding_number_yz(m, mesh, i=nx)

    v -= compute_winding_number_xz(m, mesh, j=1)
    v += compute_winding_number_xz(m, mesh, j=ny)

    return v
end


function compute_skyrmion_number(m::Array{T, 1}, mesh::Mesh) where {T<:AbstractFloat}
    nx,ny,nz = mesh.nx, mesh.ny, mesh.nz
    v = zeros(T, 3*ny*ny*nz)
    compute_skyrmion_number(v, m, mesh)
    return sum(v)
end

#compute the guiding centre, Dynamics of magnetic vortices,
#N.Papanicolaou, T.N. Tomaras 360, 425-462, (1991)
function compute_guiding_centre(m::Array{T, 1}, mesh::Mesh) where {T<:AbstractFloat}
    nx,ny,nz = mesh.nx, mesh.ny, mesh.nz
    dx, dy = mesh.dx, mesh.dy
    Rxs = zeros(nz)
    Rys = zeros(nz)
    for k = 1:nz
        sum, Rx, Ry = 0.0, 0.0, 0.0
        for j = 1:ny, i=1:nx
            id = index(i, j, k, nx, ny, nz)
            mx,my,mz = m[3*id-2],m[3*id-1],m[3*id]
            sx1,sy1,sz1 = T(0),T(0),T(0)
            sx2,sy2,sz2 = T(0),T(0),T(0)
            id1 = 3*_x_minus_one(i, id, nx, ny, nz, mesh.xperiodic)
            id2 = 3*_y_minus_one(j, id, nx, ny, nz, mesh.yperiodic)
            charge = 0
            if id1>0 && id2>0
                sx1,sy1,sz1 = m[id1-2],m[id1-1],m[id1]
                sx2,sy2,sz2 = m[id2-2],m[id2-1],m[id2]
                charge += Berg_Omega(sx2, sy2, sz2, mx, my, mz, sx1, sy1, sz1)
            end

            id1 = 3*_x_plus_one(i, id, nx, ny, nz, mesh.xperiodic)
            id2 = 3*_y_plus_one(j, id, nx, ny, nz, mesh.yperiodic)
            if id1>0 && id2>0
                sx1,sy1,sz1 = m[id1-2],m[id1-1],m[id1]
                sx2,sy2,sz2 = m[id2-2],m[id2-1],m[id2]
                charge += Berg_Omega(sx2, sy2, sz2, mx, my, mz, sx1, sy1, sz1)
            end

            sum += charge
            Rx += i * dx * charge;
            Ry += j * dy * charge;
        end
        if sum == 0.0
            sum = 1.0
        end
        Rxs[k] = Rx/sum
        Rys[k] = Ry/sum
    end
    return Rxs, Rys
end
