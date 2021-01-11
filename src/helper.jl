function  init_scalar!(v::Array{T, 1}, mesh::Mesh, init::Number) where {T<:AbstractFloat}
    for i = 1:mesh.nxyz
        v[i] = init
    end
    return nothing
end

function  init_scalar!(v::Array{T1, 1}, mesh::Mesh, init::Array{T2, 1}) where {T1,T2<:AbstractFloat}
    v .= init
    return nothing
end

function  init_scalar!(v::Array{Bool, 1}, mesh::Mesh, init::Array{Bool, 1})
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

function init_scalar!(v::Array{Bool, 1}, mesh::Mesh, init_fun::Function)
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

function init_vector!(v::Array{T, 1}, mesh::Mesh, init::Tuple{Real,Real,Real}) where {T<:AbstractFloat}
  nxyz = mesh.nx*mesh.ny*mesh.nz
  b = reshape(v, 3, nxyz)
  b[1, :] .= init[1]
  b[2, :] .= init[2]
  b[3, :] .= init[3]
  return nothing
end

function init_vector!(v::Array{T1, 1}, mesh::Mesh, init::Array{T2, 1}) where {T1,T2<:AbstractFloat}
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

function compute_dm!(dm::Array{Float64, 1}, m1::Array{Float64, 1}, m2::Array{Float64, 1}, N::Int64)
  for i = 1:N
    j = 3*i - 2
    mx = m1[j] - m2[j]
    my = m1[j+1] - m2[j+1]
    mz = m1[j+2] - m2[j+2]
    dm[i] = sqrt(mx*mx + my*my + mz*mz)
  end
  return nothing
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

#shape factor is defined as (1/4*pi) \int \partial_i m * \partial_j m dx dy
function compute_shape_factor(m::Array{T, 1}, mesh::Mesh) where {T<:AbstractFloat}
    nx,ny,nz = mesh.nx, mesh.ny, mesh.nz
    pxm, pym = partial_xy(m, mesh)
    eta_xx = 0
    eta_xy = 0
    eta_yx = 0
    eta_yy = 0
    for k = 1:nz, j = 1:ny, i=1:nx
        id = index(i, j, k, nx, ny, nz)
        px_mx,px_my,px_mz = pxm[3*id-2],pxm[3*id-1],pxm[3*id]
        py_mx,py_my,py_mz = pym[3*id-2],pym[3*id-1],pym[3*id]
        eta_xx += dot_product(px_mx,px_my,px_mz, px_mx,px_my,px_mz)
        eta_xy += dot_product(px_mx,px_my,px_mz, py_mx,py_my,py_mz)
        eta_yx += dot_product(py_mx,py_my,py_mz, px_mx,px_my,px_mz)
        eta_yy += dot_product(py_mx,py_my,py_mz, py_mx,py_my,py_mz)
    end
    factor = mesh.dx*mesh.dy/(4*pi)
    return eta_xx*factor, eta_xy*factor, eta_yx*factor, eta_yy*factor
end


function compute_winding_number_yz(v::Array{T, 1}, m::Array{T, 1}, mesh::Mesh) where {T<:AbstractFloat}
    nx,ny,nz = mesh.nx, mesh.ny, mesh.nz
    for k = 1:nz, j = 1:ny, i=1:nx
        id = index(i, j, k, nx, ny, nz)
        mx,my,mz = m[3*id-2],m[3*id-1],m[3*id]
        sx1,sy1,sz1 = T(0),T(0),T(0)
        sx2,sy2,sz2 = T(0),T(0),T(0)
        id1 = 3*_y_minus_one(j, id, nx, ny, nz, mesh.yperiodic)
        id2 = 3*_z_minus_one(k, id, nx, ny, nz, mesh.zperiodic)
        v[id] = 0
        if id1>0 && id2>0
            sx1,sy1,sz1 = m[id1-2],m[id1-1],m[id1]
            sx2,sy2,sz2 = m[id2-2],m[id2-1],m[id2]
            v[id]  += Berg_Omega(sx2, sy2, sz2, mx, my, mz, sx1, sy1, sz1)
        end

        id1 = 3*_y_plus_one(j, id, nx, ny, nz, mesh.yperiodic)
        id2 = 3*_z_plus_one(k, id, nx, ny, nz, mesh.zperiodic)
        if id1>0 && id2>0
            sx1,sy1,sz1 = m[id1-2],m[id1-1],m[id1]
            sx2,sy2,sz2 = m[id2-2],m[id2-1],m[id2]
            v[id] += Berg_Omega(sx2, sy2, sz2, mx, my, mz, sx1, sy1, sz1)
        end
        v[id] /= (4*pi);
    end

    return nothing;
end

function compute_winding_number_zx(v::Array{T, 1}, m::Array{T, 1}, mesh::Mesh) where {T<:AbstractFloat}
    nx,ny,nz = mesh.nx, mesh.ny, mesh.nz
    for k = 1:nz, j = 1:ny, i=1:nx
        id = index(i, j, k, nx, ny, nz)
        mx,my,mz = m[3*id-2],m[3*id-1],m[3*id]
        sx1,sy1,sz1 = T(0),T(0),T(0)
        sx2,sy2,sz2 = T(0),T(0),T(0)
        id1 = 3*_z_minus_one(k, id, nx, ny, nz, mesh.zperiodic)
        id2 = 3*_x_minus_one(i, id, nx, ny, nz, mesh.xperiodic)
        v[id] = 0
        if id1>0 && id2>0
            sx1,sy1,sz1 = m[id1-2],m[id1-1],m[id1]
            sx2,sy2,sz2 = m[id2-2],m[id2-1],m[id2]
            v[id] += Berg_Omega(sx2, sy2, sz2, mx, my, mz, sx1, sy1, sz1)
        end

        id1 = 3*_z_plus_one(k, id, nx, ny, nz, mesh.zperiodic)
        id2 = 3*_x_plus_one(i, id, nx, ny, nz, mesh.xperiodic)
        if id1>0 && id2>0
            sx1,sy1,sz1 = m[id1-2],m[id1-1],m[id1]
            sx2,sy2,sz2 = m[id2-2],m[id2-1],m[id2]
            v[id] += Berg_Omega(sx2, sy2, sz2, mx, my, mz, sx1, sy1, sz1)
        end
        v[id] /= (4*pi);

    end
    return nothing;
end

#We define the winding number density as
# \rho = \nabla G = \partial_x G_x + \partial_y G_y + \partial_z G_z
# where G_x, G_y, G_z are the skyrmion number defined in yz, zx, xy-plane.
function compute_winding_number_3d(v::Array{T, 1}, m::Array{T, 1}, mesh::Mesh) where {T<:AbstractFloat}
    nx,ny,nz = mesh.nx, mesh.ny, mesh.nz
    dx, dy, dz = mesh.dx, mesh.dy, mesh.dz
    ngbs = mesh.ngbs
    vx = zeros(T, ny*ny*nz)
    vy = zeros(T, ny*ny*nz)
    vz = zeros(T, ny*ny*nz)
    compute_winding_number_yz(vx, m, mesh)
    compute_winding_number_zx(vy, m, mesh)
    compute_skyrmion_number(vz, m, mesh) #compute_winding_number_xy

    for k = 1:nz, j = 1:ny, i=1:nx
        id = index(i, j, k, nx, ny, nz)
        v[id] = 0

        #x-direction
        i1 = ngbs[1,id]
        i2 = ngbs[2,id]
        factor = i1*i2>0 ? 1/(2*dx) : 1/dx
        i1 < 0 && (i1 = id)
        i2 < 0 && (i2 = id)
        v[id] += (vx[i2] - vx[i1]) * factor*dx;

        #y-direction
        i1 = ngbs[3,id]
        i2 = ngbs[4,id]
        factor = i1*i2>0 ? 1/(2*dy) : 1/dy
        i1 < 0 && (i1 = id)
        i2 < 0 && (i2 = id)
        v[id] += (vy[i2] - vy[i1]) * factor*dy;

        #z-direction
        i1 = ngbs[5,id]
        i2 = ngbs[6,id]
        factor = i1*i2>0 ? 1/(2*dz) : 1/dz
        i1 < 0 && (i1 = id)
        i2 < 0 && (i2 = id)
        v[id] += (vz[i2] - vz[i1]) * factor*dz;

    end
    return nothing
end

function winding_number_3d(m::Array{T, 1}, mesh::Mesh) where {T<:AbstractFloat}
    nx,ny,nz = mesh.nx, mesh.ny, mesh.nz
    v = zeros(T, nx*ny*nz)
    compute_winding_number_3d(v, m, mesh)
    return sum(v)
end


function compute_skyrmion_number(m::Array{T, 1}, mesh::Mesh) where {T<:AbstractFloat}
    nx,ny,nz = mesh.nx, mesh.ny, mesh.nz
    v = zeros(T, nx*ny*nz)
    compute_skyrmion_number(v, m, mesh)
    return sum(v)
end

"""
  compute_skyrmion_number_layers(ovf_name)

compute the skyrmion number of each layer of the given ovf and return an array.

  ```julia
      skx_number = compute_skyrmion_number_layers("my.ovf")
  ```
"""
function compute_skyrmion_number_layers(fname::String)
    ovf = read_ovf(fname)
    m = ovf.data
    nx,ny,nz = ovf.xnodes, ovf.ynodes, ovf.znodes
    dx,dy,dz = ovf.xstepsize, ovf.ystepsize, ovf.zstepsize

    mesh = FDMesh(nx=nx, ny=ny, nz=nz, dx=dx, dy=dy, dz=dz)

    v = zeros(nx*ny*nz)
    compute_skyrmion_number(v, m, mesh)

    b = reshape(v,(nx,ny,nz))
    skx_number = zeros(nz)
    for k =1:nz
        skx_number[k] = sum(b[:,:,k])
    end

    return skx_number
end

function compute_skyrmion_number(fname::String)
    ovf = read_ovf(fname)
    m = ovf.data
    nx,ny,nz = ovf.xnodes, ovf.ynodes, ovf.znodes
    dx,dy,dz = ovf.xstepsize, ovf.ystepsize, ovf.zstepsize

    mesh = FDMesh(nx=nx, ny=ny, nz=nz, dx=dx, dy=dy, dz=dz)

    return compute_skyrmion_number(m, mesh)
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

#F_i = \vec{p} \cdot (\vec{m} \times \partial_i \vec{m})
function compute_cpp_force(m::Array{T, 1}, mesh::Mesh; p=(0,1,0)) where {T<:AbstractFloat}
    nx,ny,nz = mesh.nx, mesh.ny, mesh.nz
    dx, dy = mesh.dx, mesh.dy
    nxyz = mesh.nxyz
    ngbs = mesh.ngbs

    pxm, pym = partial_xy(m, mesh)
    fx = zeros(T, nxyz)
    fy = zeros(T, nxyz)
    for i = 1:nxyz
      j = 3*i-2
      mx, my, mz = m[j], m[j+1], m[j+2]
      fx[i] = p[1]* cross_x(mx, my, mz, pxm[j], pxm[j+1], pxm[j+2])
      fx[i] += p[2]* cross_y(mx, my, mz, pxm[j], pxm[j+1], pxm[j+2])
      fx[i] += p[3]* cross_z(mx, my, mz, pxm[j], pxm[j+1], pxm[j+2])

      fy[i] = p[1]* cross_x(mx, my, mz, pym[j], pym[j+1], pym[j+2])
      fy[i] += p[2]* cross_y(mx, my, mz, pym[j], pym[j+1], pym[j+2])
      fy[i] += p[3]* cross_z(mx, my, mz, pym[j], pym[j+1], pym[j+2])
    end
    return fx, fy
end
