export compute_skyrmion_number, compute_guiding_center

function compute_skyrmion_number(v::Array{T,1}, m::Array{T,1},
                                 mesh::Mesh) where {T<:AbstractFloat}
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
    for k in 1:nz, j in 1:ny, i in 1:nx
        id = index(i, j, k, nx, ny, nz)
        mx, my, mz = m[3 * id - 2], m[3 * id - 1], m[3 * id]
        sx1, sy1, sz1 = T(0), T(0), T(0)
        sx2, sy2, sz2 = T(0), T(0), T(0)
        id1 = 3 * _x_minus_one(i, id, nx, ny, nz, mesh.xperiodic)
        id2 = 3 * _y_minus_one(j, id, nx, ny, nz, mesh.yperiodic)
        v[id] = 0
        if id1 > 0 && id2 > 0
            sx1, sy1, sz1 = m[id1 - 2], m[id1 - 1], m[id1]
            sx2, sy2, sz2 = m[id2 - 2], m[id2 - 1], m[id2]
            v[id] += Berg_Omega(sx2, sy2, sz2, mx, my, mz, sx1, sy1, sz1)
        end

        id1 = 3 * _x_plus_one(i, id, nx, ny, nz, mesh.xperiodic)
        id2 = 3 * _y_plus_one(j, id, nx, ny, nz, mesh.yperiodic)
        if id1 > 0 && id2 > 0
            sx1, sy1, sz1 = m[id1 - 2], m[id1 - 1], m[id1]
            sx2, sy2, sz2 = m[id2 - 2], m[id2 - 1], m[id2]
            v[id] += Berg_Omega(sx2, sy2, sz2, mx, my, mz, sx1, sy1, sz1)
        end
        v[id] /= (4 * pi)
    end
    return nothing
end

function compute_skyrmion_number(v::Array{T,1}, m::Array{T,1}, mesh::Mesh,
                                 shape::CSGShape) where {T<:AbstractFloat}
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
    dx, dy, dz = mesh.dx, mesh.dy, mesh.dz
    for k in 1:nz, j in 1:ny, i in 1:nx
        id = index(i, j, k, nx, ny, nz)
        mx, my, mz = m[3 * id - 2], m[3 * id - 1], m[3 * id]
        sx1, sy1, sz1 = T(0), T(0), T(0)
        sx2, sy2, sz2 = T(0), T(0), T(0)
        id1 = 3 * _x_minus_one(i, id, nx, ny, nz, mesh.xperiodic)
        id2 = 3 * _y_minus_one(j, id, nx, ny, nz, mesh.yperiodic)
        v[id] = 0
        x = (i - 0.5 - nx / 2) * dx
        y = (j - 0.5 - ny / 2) * dy
        z = (k - 0.5 - nz / 2) * dz
        if id1 > 0 && id2 > 0 && inside(shape, (x, y, z))
            sx1, sy1, sz1 = m[id1 - 2], m[id1 - 1], m[id1]
            sx2, sy2, sz2 = m[id2 - 2], m[id2 - 1], m[id2]
            v[id] += Berg_Omega(sx2, sy2, sz2, mx, my, mz, sx1, sy1, sz1)
        end

        id1 = 3 * _x_plus_one(i, id, nx, ny, nz, mesh.xperiodic)
        id2 = 3 * _y_plus_one(j, id, nx, ny, nz, mesh.yperiodic)
        if id1 > 0 && id2 > 0 && inside(shape, (x, y, z))
            sx1, sy1, sz1 = m[id1 - 2], m[id1 - 1], m[id1]
            sx2, sy2, sz2 = m[id2 - 2], m[id2 - 1], m[id2]
            v[id] += Berg_Omega(sx2, sy2, sz2, mx, my, mz, sx1, sy1, sz1)
        end
        v[id] /= (4 * pi)
    end
    return nothing
end

function compute_skyrmion_number(m::Array{T,1}, mesh::Mesh) where {T<:AbstractFloat}
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
    v = zeros(T, nx * ny * nz)
    compute_skyrmion_number(v, m, mesh)
    return sum(v)
end

function compute_skyrmion_number(m::Array{T,1}, mesh::Mesh,
                                 shape::CSGShape) where {T<:AbstractFloat}
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
    v = zeros(T, nx * ny * nz)
    compute_skyrmion_number(v, m, mesh, shape)
    return sum(v)
end

"""
  compute_skyrmion_number_layers(fname::String)

compute the skyrmion number of each layer of the given ovf and return an array.

  ```julia
      skx_number = compute_skyrmion_number_layers("my.ovf")
  ```
"""
function compute_skyrmion_number_layers(fname::String)
    ovf = read_ovf(fname)
    m = ovf.data
    nx, ny, nz = ovf.xnodes, ovf.ynodes, ovf.znodes
    dx, dy, dz = ovf.xstepsize, ovf.ystepsize, ovf.zstepsize

    mesh = FDMesh(; nx=nx, ny=ny, nz=nz, dx=dx, dy=dy, dz=dz)

    v = zeros(nx * ny * nz)
    compute_skyrmion_number(v, m, mesh)

    b = reshape(v, (nx, ny, nz))
    skx_number = zeros(nz)
    for k in 1:nz
        skx_number[k] = sum(b[:, :, k])
    end

    return skx_number
end
"""
  compute_skyrmion_number(fname::String)

compute the total skyrmion number of the given ovf file and return a number.

  ```julia
      skx_number = compute_skyrmion_number("my.ovf")
  ```
"""
function compute_skyrmion_number(fname::String)
    ovf = read_ovf(fname)
    m = ovf.data
    nx, ny, nz = ovf.xnodes, ovf.ynodes, ovf.znodes
    dx, dy, dz = ovf.xstepsize, ovf.ystepsize, ovf.zstepsize

    mesh = FDMesh(; nx=nx, ny=ny, nz=nz, dx=dx, dy=dy, dz=dz)

    return compute_skyrmion_number(m, mesh)
end

"""
    compute_guiding_center(m::Array{T, 1}, mesh::Mesh; xmin=1, xmax=-1, ymin = 1, ymax=-1, z=1)

compute the guiding center, see [1].

[1] Dynamics of magnetic vortices, N.Papanicolaou, T.N. Tomaras 360, 425-462, (1991)
"""
function compute_guiding_center(m::Array{T,1}, mesh::Mesh; xmin=1, xmax=-1, ymin=1, ymax=-1,
                                z=1) where {T<:AbstractFloat}
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
    dx, dy = mesh.dx, mesh.dy

    xmax < 0 && (xmax = nx)
    ymax < 0 && (ymax = ny)

    total_charge, Rx, Ry = 0.0, 0.0, 0.0
    for j in ymin:ymax, i in xmin:xmax
        id = index(i, j, z, nx, ny, nz)
        mx, my, mz = m[3 * id - 2], m[3 * id - 1], m[3 * id]
        sx1, sy1, sz1 = T(0), T(0), T(0)
        sx2, sy2, sz2 = T(0), T(0), T(0)
        id1 = 3 * _x_minus_one(i, id, nx, ny, nz, mesh.xperiodic)
        id2 = 3 * _y_minus_one(j, id, nx, ny, nz, mesh.yperiodic)
        charge = 0
        if id1 > 0 && id2 > 0
            sx1, sy1, sz1 = m[id1 - 2], m[id1 - 1], m[id1]
            sx2, sy2, sz2 = m[id2 - 2], m[id2 - 1], m[id2]
            charge += Berg_Omega(sx2, sy2, sz2, mx, my, mz, sx1, sy1, sz1)
        end

        id1 = 3 * _x_plus_one(i, id, nx, ny, nz, mesh.xperiodic)
        id2 = 3 * _y_plus_one(j, id, nx, ny, nz, mesh.yperiodic)
        if id1 > 0 && id2 > 0
            sx1, sy1, sz1 = m[id1 - 2], m[id1 - 1], m[id1]
            sx2, sy2, sz2 = m[id2 - 2], m[id2 - 1], m[id2]
            charge += Berg_Omega(sx2, sy2, sz2, mx, my, mz, sx1, sy1, sz1)
        end

        total_charge += charge
        Rx += i * dx * charge
        Ry += j * dy * charge
    end

    if total_charge == 0.0
        total_charge = 1.0
    end
    return Rx / total_charge, Ry / total_charge
end

"""
    function compute_guiding_center(sim::AbstractSim; xmin=1, xmax=-1, ymin = 1, ymax=-1, z=1)

compute the guiding center.
"""
function compute_guiding_center(sim::AbstractSim; xmin=1, xmax=-1, ymin=1, ymax=-1, z=1)
    spin = Array(sim.spin)
    mesh = sim.mesh

    if isa(mesh, CylindricalTubeMesh)
        return compute_guiding_center(sim, sim.mesh)
    end

    Rx, Ry = compute_guiding_center(spin, mesh; xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
                                    z=z)
    return Rx, Ry
end

#shape factor is defined as (1/4*pi) \int \partial_i m * \partial_j m dx dy
function compute_shape_factor(m::Array{T,1}, mesh::Mesh; r=0) where {T<:AbstractFloat}
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
    pxm, pym = partial_xy(m, mesh)
    eta_xx = 0
    eta_xy = 0
    eta_yx = 0
    eta_yy = 0
    if r == 0
        r = max(nx, ny)
    end
    x_c = nx / 2.0 + 0.5
    y_c = ny / 2.0 + 0.5
    factor = mesh.dx * mesh.dy / (4 * pi * nz)
    for k in 1:nz, j in 1:ny, i in 1:nx
        id = index(i, j, k, nx, ny, nz)
        px_mx, px_my, px_mz = pxm[3 * id - 2], pxm[3 * id - 1], pxm[3 * id]
        py_mx, py_my, py_mz = pym[3 * id - 2], pym[3 * id - 1], pym[3 * id]
        if (i - x_c)^2 + (j - y_c)^2 < r^2
            eta_xx += dot_product(px_mx, px_my, px_mz, px_mx, px_my, px_mz)
            eta_xy += dot_product(px_mx, px_my, px_mz, py_mx, py_my, py_mz)
            eta_yx += dot_product(py_mx, py_my, py_mz, px_mx, px_my, px_mz)
            eta_yy += dot_product(py_mx, py_my, py_mz, py_mx, py_my, py_mz)
        end
    end

    return eta_xx * factor, eta_xy * factor, eta_yx * factor, eta_yy * factor
end

function compute_shape_factor(m::Array{T,1}, x_c::Int64, y_c::Int64, mesh::Mesh;
                              r=0) where {T<:AbstractFloat}
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
    pxm, pym = partial_xy(m, mesh)
    eta_xx = 0
    eta_xy = 0
    eta_yx = 0
    eta_yy = 0
    if r == 0
        r = max(nx, ny)
    end
    #x_c = nx/2.0 + 0.5
    #y_c = ny/2.0 + 0.5
    factor = mesh.dx * mesh.dy / (4 * pi * nz)
    for k in 1:nz, j in 1:ny, i in 1:nx
        id = index(i, j, k, nx, ny, nz)
        px_mx, px_my, px_mz = pxm[3 * id - 2], pxm[3 * id - 1], pxm[3 * id]
        py_mx, py_my, py_mz = pym[3 * id - 2], pym[3 * id - 1], pym[3 * id]
        if (i - x_c)^2 + (j - y_c)^2 < r^2
            eta_xx += dot_product(px_mx, px_my, px_mz, px_mx, px_my, px_mz)
            eta_xy += dot_product(px_mx, px_my, px_mz, py_mx, py_my, py_mz)
            eta_yx += dot_product(py_mx, py_my, py_mz, px_mx, px_my, px_mz)
            eta_yy += dot_product(py_mx, py_my, py_mz, py_mx, py_my, py_mz)
        end
    end

    return eta_xx * factor, eta_xy * factor, eta_yx * factor, eta_yy * factor
end

# The tensor B is defined as
# \mathcal{B}_{i j}=\int\left(\mathbf{m} \times \partial_i \mathbf{m}\right)_j d x d y
function compute_tensor_B(m::Array{T,1}, mesh::Mesh; r=-1) where {T<:AbstractFloat}
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
    pxm, pym = partial_xy(m, mesh)
    eta_xx = 0
    eta_xy = 0
    eta_yx = 0
    eta_yy = 0
    if r < 0
        r = max(nx, ny)
    end
    x_c = nx / 2.0 + 0.5
    y_c = ny / 2.0 + 0.5
    factor = mesh.dx * mesh.dy / nz
    for k in 1:nz, j in 1:ny, i in 1:nx
        id = index(i, j, k, nx, ny, nz)

        s = 3 * id - 2
        mx, my, mz = m[s], m[s + 1], m[s + 2]

        if (i - x_c)^2 + (j - y_c)^2 < r^2
            eta_xx += cross_x(mx, my, mz, pxm[s], pxm[s + 1], pxm[s + 2])
            eta_xy += cross_y(mx, my, mz, pxm[s], pxm[s + 1], pxm[s + 2])
            eta_yx += cross_x(mx, my, mz, pym[s], pym[s + 1], pym[s + 2])
            eta_yy += cross_y(mx, my, mz, pym[s], pym[s + 1], pym[s + 2])
        end
    end

    return eta_xx * factor, eta_xy * factor, eta_yx * factor, eta_yy * factor
end
