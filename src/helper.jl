
export compute_skyrmion_number, compute_guiding_center

function init_scalar!(v::AbstractArray{T,1}, mesh::Mesh,
                      init::Number) where {T<:AbstractFloat}
    a = isa(v, Array) ? v : Array(v)
    for i in 1:(mesh.n_total)
        a[i] = init
    end
    isa(v, Array) || copyto!(v, a)
    return true
end

function init_scalar!(v::AbstractArray{T1,1}, mesh::Mesh, init::Array{T2}) where {T1,T2}
    a = isa(v, Array) ? v : Array(v)
    a .= init
    isa(v, Array) || copyto!(v, a)
    return true
end

function init_scalar!(v::AbstractArray{T,1}, mesh::Mesh, init_fun::Function) where {T}
    a = isa(v, Array) ? v : Array(v)
    for k in 1:(mesh.nz), j in 1:(mesh.ny), i in 1:(mesh.nx)
        id = index(i, j, k, mesh.nx, mesh.ny, mesh.nz)
        a[id] = init_fun(i, j, k, mesh.dx, mesh.dy, mesh.dz)
    end
    isa(v, Array) || copyto!(v, a)
    return true
end

function init_scalar!(v::AbstractArray{T,1}, mesh::Mesh, shape::Union{CSGNode, Shape},
                      init::Number) where {T}
    a = isa(v, Array) ? v : Array(v)
    dx, dy, dz = mesh.dx, mesh.dy, mesh.dz
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
    for k in 1:nz, j in 1:ny, i in 1:nx
        id = index(i, j, k, mesh.nx, mesh.ny, mesh.nz)
        x = (i - 0.5 - nx / 2) * dx
        y = (j - 0.5 - ny / 2) * dy
        z = (k - 0.5 - nz / 2) * dz
        if point_inside_shape((x, y, z), shape)
            a[id] = init
        end
    end
    isa(v, Array) || copyto!(v, a)
    return true
end

function init_vector!(v::Array{T,1}, mesh::Mesh, init::Function) where {T<:AbstractFloat}
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
    n_total = nx * ny * nz
    dx, dy, dz = mesh.dx, mesh.dy, mesh.dz
    b = reshape(v, 3, n_total)

    nargs = methods(init)[1].nargs - 1
    if nargs == 6
        for i in 1:nx, j in 1:ny, k in 1:nz
            id = index(i, j, k, nx, ny, nz)
            vec_value = init(i, j, k, dx, dy, dz)
            # ignore the values for specfic positions that the user do not want to provide or change.
            if vec_value !== nothing
                b[:, id] .= vec_value[:]
            end
        end

    elseif nargs == 3
        for i in 1:nx, j in 1:ny, k in 1:nz
            id = index(i, j, k, nx, ny, nz)
            x = (i - 0.5 - nx / 2) * dx
            y = (j - 0.5 - ny / 2) * dy
            z = (k - 0.5 - nz / 2) * dz

            vec_value = init(x, y, z)
            # ignore the values for specfic positions that the user do not want to provide or change.
            if vec_value !== nothing
                b[:, id] .= vec_value[:]
            end
        end
    else
        error("The input function should have either 6 or 3 arguments.")
    end
    
    if NaN in v
        error("NaN is given by the input function.")
    end
    return nothing
end

function init_vector!(v::Array{T,1}, mesh::Mesh,
                      init::Tuple{Real,Real,Real}) where {T<:AbstractFloat}
    #n_total = mesh.nx*mesh.ny*mesh.nz
    N = length(v)
    b = reshape(v, 3, div(N, 3))
    b[1, :] .= init[1]
    b[2, :] .= init[2]
    b[3, :] .= init[3]
    return nothing
end

function init_vector!(v::Array{T1,1}, mesh::Mesh,
                      init::Array{T2,1}) where {T1,T2<:AbstractFloat}
    v .= init
    return nothing
end

function normalise(a::AbstractArray{T,1}, N::Int64) where {T<:AbstractFloat}
    @kernel function local_kernal!(a)
        id = @index(Global)
        j = 3 * id - 2

        @inbounds m2 = a[j] * a[j] + a[j + 1] * a[j + 1] + a[j + 2] * a[j + 2]
        if m2 > 0
            length::T = 1 / sqrt(m2)
            @inbounds a[j] *= length
            @inbounds a[j + 1] *= length
            @inbounds a[j + 2] *= length
        end
    end

    local_kernal!(get_backend(a), groupsize[])(a; ndrange=N)
    return nothing
end

function compute_dm!(dm::AbstractArray{T,1}, m1::AbstractArray{T,1}, m2::AbstractArray{T,1},
                     N::Int64) where {T<:AbstractFloat}
    @kernel function local_kernal!(c, a, b)
        I = @index(Global)
        j = 3 * I - 2
        @inbounds mx = a[j] - b[j]
        @inbounds my = a[j + 1] - b[j + 1]
        @inbounds mz = a[j + 2] - b[j + 2]
        @inbounds c[I] = sqrt(mx * mx + my * my + mz * mz)
    end

    local_kernal!(get_backend(dm), groupsize[])(dm, m1, m2; ndrange=N)
    return nothing
end

function omega_to_spin(omega::AbstractArray{T,1}, spin::AbstractArray{T,1},
                       spin_next::AbstractArray{T,1}, N::Int64) where {T<:AbstractFloat}
    #compute Cay(Omega).m where Cay(Omega) = (I - 1/2 Omega)^-1 (I + 1/2 Omega)
    #where Omega = Skew[w1, w2, w3] = {{0, -w3, w2}, {w3, 0, -w1}, {-w2, w1, 0}}
    @kernel function local_kernal!(a, b, c)
        I = @index(Global)
        j = 3 * I - 2
        @inbounds w1 = a[j] * 0.5
        @inbounds w2 = a[j + 1] * 0.5
        @inbounds w3 = a[j + 2] * 0.5
        @inbounds m1 = b[j]
        @inbounds m2 = b[j + 1]
        @inbounds m3 = b[j + 2]
        r = 1 + w1 * w1 + w2 * w2 + w3 * w3
        a11 = 1 + w1 * w1 - w2 * w2 - w3 * w3
        a12 = 2 * (w1 * w2 - w3)
        a13 = 2 * (w2 + w1 * w3)
        a21 = 2 * (w1 * w2 + w3)
        a22 = 1 - w1 * w1 + w2 * w2 - w3 * w3
        a23 = -2 * (w1 - w2 * w3)
        a31 = 2 * (-w2 + w1 * w3)
        a32 = 2 * (w1 + w2 * w3)
        a33 = 1 - w1 * w1 - w2 * w2 + w3 * w3
        @inbounds c[j] = (a11 * m1 + a12 * m2 + a13 * m3) / r
        @inbounds c[j + 1] = (a21 * m1 + a22 * m2 + a23 * m3) / r
        @inbounds c[j + 2] = (a31 * m1 + a32 * m2 + a33 * m3) / r
    end

    local_kernal!(default_backend[], groupsize[])(omega, spin, spin_next; ndrange=N)
    return nothing
end

function error_length_m(a::Array{Float64,1}, N::Int64)
    maxlength = 0.0
    minlength = 1.0
    for i in 0:(N - 1)
        j = 3 * i + 1
        length = sqrt(a[j] * a[j] + a[j + 1] * a[j + 1] + a[j + 2] * a[j + 2])
        if length > maxlength
            maxlength = length
        end
        if length < minlength
            minlength = length
        end
    end
    return maxlength - minlength
end

function omega_to_spin(omega::Array{Float64,1}, spin::Array{Float64,1},
                       spin_next::Array{Float64,1}, N::Int64)
    #compute Cay(Omega).m where Cay(Omega) = (I - 1/2 Omega)^-1 (I + 1/2 Omega)
    #where Omega = Skew[w1, w2, w3] = {{0, -w3, w2}, {w3, 0, -w1}, {-w2, w1, 0}}
    for i in 0:(N - 1)
        j = 3 * i + 1
        w1 = omega[j] * 0.5
        w2 = omega[j + 1] * 0.5
        w3 = omega[j + 2] * 0.5
        m1 = spin[j]
        m2 = spin[j + 1]
        m3 = spin[j + 2]
        r = 1 + w1 * w1 + w2 * w2 + w3 * w3
        a11 = 1 + w1 * w1 - w2 * w2 - w3 * w3
        a12 = 2 * (w1 * w2 - w3)
        a13 = 2 * (w2 + w1 * w3)
        a21 = 2 * (w1 * w2 + w3)
        a22 = 1 - w1 * w1 + w2 * w2 - w3 * w3
        a23 = -2 * (w1 - w2 * w3)
        a31 = 2 * (-w2 + w1 * w3)
        a32 = 2 * (w1 + w2 * w3)
        a33 = 1 - w1 * w1 - w2 * w2 + w3 * w3
        spin_next[j] = (a11 * m1 + a12 * m2 + a13 * m3) / r
        spin_next[j + 1] = (a21 * m1 + a22 * m2 + a23 * m3) / r
        spin_next[j + 2] = (a31 * m1 + a32 * m2 + a33 * m3) / r
    end
end
function compute_error2(error::Array{Float64,1}, N::Int64)
    norm = 0.0
    for i in 1:N
        norm += error[i]^2
    end
    return sqrt(norm / N)
end

function compute_dmdt(m1::Array{Float64,1}, m2::Array{Float64,1}, N::Int64, dt::Float64)
    max_dmdt = 0.0
    for i in 0:(N - 1)
        j = 3 * i + 1
        mx = m1[j] - m2[j]
        my = m1[j + 1] - m2[j + 1]
        mz = m1[j + 2] - m2[j + 2]
        dmdt = sqrt(mx * mx + my * my + mz * mz)
        if dmdt > max_dmdt
            max_dmdt = dmdt
        end
    end
    return max_dmdt / dt
end

function compute_dm!(dm::Array{Float64,1}, m1::Array{Float64,1}, m2::Array{Float64,1},
                     N::Int64)
    for i in 1:N
        j = 3 * i - 2
        mx = m1[j] - m2[j]
        my = m1[j + 1] - m2[j + 1]
        mz = m1[j + 2] - m2[j + 2]
        dm[i] = sqrt(mx * mx + my * my + mz * mz)
    end
    return nothing
end

function compute_dm_step(m1::Array{Float64,1}, m2::Array{Float64,1}, N::Int64)
    max_dm = 0.0
    for i in 0:(N - 1)
        j = 3 * i + 1
        mx = m1[j] - m2[j]
        my = m1[j + 1] - m2[j + 1]
        mz = m1[j + 2] - m2[j + 2]
        dmdt = sqrt(mx * mx + my * my + mz * mz)
        if dmdt > max_dm
            max_dm = dmdt
        end
    end
    return max_dm
end

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

function compute_skyrmion_number(v::Array{T,1}, m::Array{T,1}, mesh::Mesh, shape::Union{CSGNode, Shape}) where {T<:AbstractFloat}
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
        if id1 > 0 && id2 > 0 && point_inside_shape((x, y, z), shape)
            sx1, sy1, sz1 = m[id1 - 2], m[id1 - 1], m[id1]
            sx2, sy2, sz2 = m[id2 - 2], m[id2 - 1], m[id2]
            v[id] += Berg_Omega(sx2, sy2, sz2, mx, my, mz, sx1, sy1, sz1)
        end

        id1 = 3 * _x_plus_one(i, id, nx, ny, nz, mesh.xperiodic)
        id2 = 3 * _y_plus_one(j, id, nx, ny, nz, mesh.yperiodic)
        if id1 > 0 && id2 > 0 && point_inside_shape((x, y, z), shape)
            sx1, sy1, sz1 = m[id1 - 2], m[id1 - 1], m[id1]
            sx2, sy2, sz2 = m[id2 - 2], m[id2 - 1], m[id2]
            v[id] += Berg_Omega(sx2, sy2, sz2, mx, my, mz, sx1, sy1, sz1)
        end
        v[id] /= (4 * pi)
    end
    return nothing
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

function compute_winding_number_yz(v::Array{T,1}, m::Array{T,1},
                                   mesh::Mesh) where {T<:AbstractFloat}
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
    for k in 1:nz, j in 1:ny, i in 1:nx
        id = index(i, j, k, nx, ny, nz)
        mx, my, mz = m[3 * id - 2], m[3 * id - 1], m[3 * id]
        sx1, sy1, sz1 = T(0), T(0), T(0)
        sx2, sy2, sz2 = T(0), T(0), T(0)
        id1 = 3 * _y_minus_one(j, id, nx, ny, nz, mesh.yperiodic)
        id2 = 3 * _z_minus_one(k, id, nx, ny, nz, mesh.zperiodic)
        v[id] = 0
        if id1 > 0 && id2 > 0
            sx1, sy1, sz1 = m[id1 - 2], m[id1 - 1], m[id1]
            sx2, sy2, sz2 = m[id2 - 2], m[id2 - 1], m[id2]
            v[id] += Berg_Omega(sx2, sy2, sz2, mx, my, mz, sx1, sy1, sz1)
        end

        id1 = 3 * _y_plus_one(j, id, nx, ny, nz, mesh.yperiodic)
        id2 = 3 * _z_plus_one(k, id, nx, ny, nz, mesh.zperiodic)
        if id1 > 0 && id2 > 0
            sx1, sy1, sz1 = m[id1 - 2], m[id1 - 1], m[id1]
            sx2, sy2, sz2 = m[id2 - 2], m[id2 - 1], m[id2]
            v[id] += Berg_Omega(sx2, sy2, sz2, mx, my, mz, sx1, sy1, sz1)
        end
        v[id] /= (4 * pi)
    end

    return nothing
end

function compute_winding_number_zx(v::Array{T,1}, m::Array{T,1},
                                   mesh::Mesh) where {T<:AbstractFloat}
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
    for k in 1:nz, j in 1:ny, i in 1:nx
        id = index(i, j, k, nx, ny, nz)
        mx, my, mz = m[3 * id - 2], m[3 * id - 1], m[3 * id]
        sx1, sy1, sz1 = T(0), T(0), T(0)
        sx2, sy2, sz2 = T(0), T(0), T(0)
        id1 = 3 * _z_minus_one(k, id, nx, ny, nz, mesh.zperiodic)
        id2 = 3 * _x_minus_one(i, id, nx, ny, nz, mesh.xperiodic)
        v[id] = 0
        if id1 > 0 && id2 > 0
            sx1, sy1, sz1 = m[id1 - 2], m[id1 - 1], m[id1]
            sx2, sy2, sz2 = m[id2 - 2], m[id2 - 1], m[id2]
            v[id] += Berg_Omega(sx2, sy2, sz2, mx, my, mz, sx1, sy1, sz1)
        end

        id1 = 3 * _z_plus_one(k, id, nx, ny, nz, mesh.zperiodic)
        id2 = 3 * _x_plus_one(i, id, nx, ny, nz, mesh.xperiodic)
        if id1 > 0 && id2 > 0
            sx1, sy1, sz1 = m[id1 - 2], m[id1 - 1], m[id1]
            sx2, sy2, sz2 = m[id2 - 2], m[id2 - 1], m[id2]
            v[id] += Berg_Omega(sx2, sy2, sz2, mx, my, mz, sx1, sy1, sz1)
        end
        v[id] /= (4 * pi)
    end
    return nothing
end

#We define the winding number density as
# \rho = \nabla G = \partial_x G_x + \partial_y G_y + \partial_z G_z
# where G_x, G_y, G_z are the skyrmion number defined in yz, zx, xy-plane.
function compute_winding_number_3d(v::Array{T,1}, m::Array{T,1},
                                   mesh::Mesh) where {T<:AbstractFloat}
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
    dx, dy, dz = mesh.dx, mesh.dy, mesh.dz
    ngbs = mesh.ngbs
    vx = zeros(T, ny * ny * nz)
    vy = zeros(T, ny * ny * nz)
    vz = zeros(T, ny * ny * nz)
    compute_winding_number_yz(vx, m, mesh)
    compute_winding_number_zx(vy, m, mesh)
    compute_skyrmion_number(vz, m, mesh) #compute_winding_number_xy

    for k in 1:nz, j in 1:ny, i in 1:nx
        id = index(i, j, k, nx, ny, nz)
        v[id] = 0

        #x-direction
        i1 = ngbs[1, id]
        i2 = ngbs[2, id]
        factor = i1 * i2 > 0 ? 1 / (2 * dx) : 1 / dx
        i1 < 0 && (i1 = id)
        i2 < 0 && (i2 = id)
        v[id] += (vx[i2] - vx[i1]) * factor * dx

        #y-direction
        i1 = ngbs[3, id]
        i2 = ngbs[4, id]
        factor = i1 * i2 > 0 ? 1 / (2 * dy) : 1 / dy
        i1 < 0 && (i1 = id)
        i2 < 0 && (i2 = id)
        v[id] += (vy[i2] - vy[i1]) * factor * dy

        #z-direction
        i1 = ngbs[5, id]
        i2 = ngbs[6, id]
        factor = i1 * i2 > 0 ? 1 / (2 * dz) : 1 / dz
        i1 < 0 && (i1 = id)
        i2 < 0 && (i2 = id)
        v[id] += (vz[i2] - vz[i1]) * factor * dz
    end
    return nothing
end

function winding_number_3d(m::Array{T,1}, mesh::Mesh) where {T<:AbstractFloat}
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
    v = zeros(T, nx * ny * nz)
    compute_winding_number_3d(v, m, mesh)
    return sum(v)
end

function compute_skyrmion_number(m::Array{T,1}, mesh::Mesh) where {T<:AbstractFloat}
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
    v = zeros(T, nx * ny * nz)
    compute_skyrmion_number(v, m, mesh)
    return sum(v)
end

function compute_skyrmion_number(m::Array{T,1}, mesh::Mesh, shape::Union{CSGNode, Shape}) where {T<:AbstractFloat}
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

#F_i = \vec{p} \cdot (\vec{m} \times \partial_i \vec{m})
function compute_cpp_force(m::Array{T,1}, mesh::Mesh; p=(0, 1, 0)) where {T<:AbstractFloat}
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
    dx, dy = mesh.dx, mesh.dy
    n_total = mesh.n_total
    ngbs = mesh.ngbs

    pxm, pym = partial_xy(m, mesh)
    fx = zeros(T, n_total)
    fy = zeros(T, n_total)
    for i in 1:n_total
        j = 3 * i - 2
        mx, my, mz = m[j], m[j + 1], m[j + 2]
        fx[i] = p[1] * cross_x(mx, my, mz, pxm[j], pxm[j + 1], pxm[j + 2])
        fx[i] += p[2] * cross_y(mx, my, mz, pxm[j], pxm[j + 1], pxm[j + 2])
        fx[i] += p[3] * cross_z(mx, my, mz, pxm[j], pxm[j + 1], pxm[j + 2])

        fy[i] = p[1] * cross_x(mx, my, mz, pym[j], pym[j + 1], pym[j + 2])
        fy[i] += p[2] * cross_y(mx, my, mz, pym[j], pym[j + 1], pym[j + 2])
        fy[i] += p[3] * cross_z(mx, my, mz, pym[j], pym[j + 1], pym[j + 2])
    end
    return fx, fy
end

function extract_sweep_keys(dict::Dict)
    result = Dict{Symbol,Any}()

    for (key, value) in dict
        if endswith(String(key), "_sweep")
            new_key = Symbol(replace(String(key), "_sweep" => ""))
            result[new_key] = value
            delete!(dict, key)
        end

        if endswith(String(key), "_s") && String(key) != "mu_s"
            new_key = Symbol(replace(String(key), "_s" => ""))
            result[new_key] = value
            delete!(dict, key)
        end
    end

    return result
end

function check_sweep_lengths(dict::Dict)
    range_keys = filter(key -> endswith(String(key), "_sweep") ||
                            (endswith(String(key), "_s") && String(key) != "mu_s"),
                        keys(dict))
    lengths = [length(dict[key]) for key in range_keys]

    if length(lengths) > 1 && length(unique(lengths)) > 1
        throw(ErrorException("Error: Not all _range arrays have the same length."))
    end

    if length(lengths) == 0
        return 0
    end
    return lengths[1]
end
