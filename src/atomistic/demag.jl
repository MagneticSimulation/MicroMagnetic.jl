using FFTW
using LinearAlgebra

function init_demag(sim::AtomisticSim, Nx::Int, Ny::Int, Nz::Int)
    mesh = sim.mesh
    #max_size = max(mesh.dx, mesh.dy, mesh.dz)
    dx = mesh.dx * 1e9
    dy = mesh.dy * 1e9
    dz = mesh.dz * 1e9

    nx = mesh.nx
    ny = mesh.ny
    nz = mesh.nz

    cn = 3
    nx_fft = mesh.nx > cn ? 2 * mesh.nx : 2 * mesh.nx - 1
    ny_fft = mesh.ny > cn ? 2 * mesh.ny : 2 * mesh.ny - 1
    nz_fft = mesh.nz > cn ? 2 * mesh.nz : 2 * mesh.nz - 1

    mx_pad = create_zeros(nx_fft, ny_fft, nz_fft)
    my_pad = create_zeros(nx_fft, ny_fft, nz_fft)
    mz_pad = create_zeros(nx_fft, ny_fft, nz_fft)
    plan = plan_rfft(mx_pad)

    tensor = create_zeros(nx, ny, nz)

    #Nxx
    compute_dipolar_tensors(tensor, dipolar_tensors_kernel_xx!, Nx, Ny, Nz, dx, dy, dz)
    fill_tensors(mx_pad, tensor, false, false, false)
    tensor_xx = real(plan * mx_pad)

    #Nyy
    compute_dipolar_tensors(tensor, dipolar_tensors_kernel_yy!, Nx, Ny, Nz, dx, dy, dz)
    fill_tensors(mx_pad, tensor, false, false, false)
    tensor_yy = real(plan * mx_pad)

    #Nzz
    compute_dipolar_tensors(tensor, dipolar_tensors_kernel_zz!, Nx, Ny, Nz, dx, dy, dz)
    fill_tensors(mx_pad, tensor, false, false, false)
    tensor_zz = real(plan * mx_pad)

    #Nxy
    compute_dipolar_tensors(tensor, dipolar_tensors_kernel_xy!, Nx, Ny, Nz, dx, dy, dz)
    fill_tensors(mx_pad, tensor, true, true, false)
    tensor_xy = real(plan * mx_pad)

    #Nxz
    compute_dipolar_tensors(tensor, dipolar_tensors_kernel_xz!, Nx, Ny, Nz, dx, dy, dz)
    fill_tensors(mx_pad, tensor, true, false, true)
    tensor_xz = real(plan * mx_pad)

    #Nyz
    compute_dipolar_tensors(tensor, dipolar_tensors_kernel_yz!, Nx, Ny, Nz, dx, dy, dz)
    fill_tensors(mx_pad, tensor, false, true, true)
    tensor_yz = real(plan * mx_pad)

    lenx = (nx_fft % 2 > 0) ? nx : nx + 1

    T = Float[]
    Mx = create_zeros(Complex{T}, lenx, ny_fft, nz_fft)
    My = create_zeros(Complex{T}, lenx, ny_fft, nz_fft)
    Mz = create_zeros(Complex{T}, lenx, ny_fft, nz_fft)

    #m_plan = plan_rfft(mx_pad)
    h_plan = plan_irfft(Mx, nx_fft)

    field = create_zeros(3 * sim.n_total)
    energy = create_zeros(sim.n_total)

    demag = Demag(nx_fft, ny_fft, nz_fft, tensor_xx, tensor_yy, tensor_zz, tensor_xy,
                  tensor_xz, tensor_yz, mx_pad, my_pad, mz_pad, Mx, My, Mz,
                  plan, h_plan, field, energy, "Demag")
    return demag
end

function effective_field(demag::Demag, sim::AtomisticSim, spin::AbstractArray{T,1},
                         t::Float64) where {T<:AbstractFloat}
    mesh = sim.mesh
    nx, ny, nz = mesh.nx, sim.mesh.ny, sim.mesh.nz

    fill!(demag.mx, 0)
    fill!(demag.my, 0)
    fill!(demag.mz, 0)

    distribute_m_atomistic(spin, demag.mx, demag.my, demag.mz, sim.mu_s, nx, ny, nz)

    #synchronize()
    mul!(demag.Mx, demag.m_plan, demag.mx)
    mul!(demag.My, demag.m_plan, demag.my)
    mul!(demag.Mz, demag.m_plan, demag.mz)

    add_tensor_M(demag.Mx, demag.My, demag.Mz, demag.tensor_xx, demag.tensor_yy,
                 demag.tensor_zz, demag.tensor_xy, demag.tensor_xz, demag.tensor_yz)
    #synchronize()

    mul!(demag.mx, demag.h_plan, demag.Mx)
    mul!(demag.my, demag.h_plan, demag.My)
    mul!(demag.mz, demag.h_plan, demag.Mz)

    collect_h_atomistic_energy(demag.field, demag.energy, spin, demag.mx, demag.my,
                               demag.mz, sim.mu_s, nx, ny, nz)

    return nothing
end

function dipolar_tensor_xx(x::Float64, y::Float64, z::Float64)
    R = x * x + y * y + z * z
    if R == 0
        return 0.0
    else
        return -(2 * x * x - y * y - z * z) / (R * R * sqrt(R))
    end
end

function dipolar_tensor_xy(x::Float64, y::Float64, z::Float64)
    R = x * x + y * y + z * z
    if R == 0
        return 0.0
    else
        return -3 * x * y / (R * R * sqrt(R))
    end
end

function dipolar_tensor_yy(x::Float64, y::Float64, z::Float64)
    return dipolar_tensor_xx(y, x, z)
end

function dipolar_tensor_zz(x::Float64, y::Float64, z::Float64)
    return dipolar_tensor_xx(z, y, x)
end

function dipolar_tensor_xz(x::Float64, y::Float64, z::Float64)
    return dipolar_tensor_xy(x, z, y)
end

function dipolar_tensor_yz(x::Float64, y::Float64, z::Float64)
    return dipolar_tensor_xy(y, z, x)
end

@kernel function dipolar_tensors_kernel_xx!(tensor, dx::Float64, dy::Float64, dz::Float64,
                                            Nx::Int64, Ny::Int64, Nz::Int64)
    nx, ny, nz = size(tensor)
    i, j, k = @index(Global, NTuple)

    sum = 0.0
    for p in (-Nx):Nx, q in (-Ny):Ny, s in (-Nz):Nz
        x = (i - 1 + p * nx) * dx
        y = (j - 1 + q * ny) * dy
        z = (k - 1 + s * nz) * dz
        sum += dipolar_tensor_xx(x, y, z)
    end
    @inbounds tensor[i, j, k] = sum
end

@kernel function dipolar_tensors_kernel_yy!(tensor, dx::Float64, dy::Float64, dz::Float64,
                                            Nx::Int64, Ny::Int64, Nz::Int64)
    nx, ny, nz = size(tensor)
    i, j, k = @index(Global, NTuple)

    sum = 0.0
    for p in (-Nx):Nx, q in (-Ny):Ny, s in (-Nz):Nz
        x = (i - 1 + p * nx) * dx
        y = (j - 1 + q * ny) * dy
        z = (k - 1 + s * nz) * dz
        sum += dipolar_tensor_yy(x, y, z)
    end
    @inbounds tensor[i, j, k] = sum
end

@kernel function dipolar_tensors_kernel_zz!(tensor, dx::Float64, dy::Float64, dz::Float64,
                                            Nx::Int64, Ny::Int64, Nz::Int64)
    nx, ny, nz = size(tensor)
    i, j, k = @index(Global, NTuple)

    sum = 0.0
    for p in (-Nx):Nx, q in (-Ny):Ny, s in (-Nz):Nz
        x = (i - 1 + p * nx) * dx
        y = (j - 1 + q * ny) * dy
        z = (k - 1 + s * nz) * dz
        sum += dipolar_tensor_zz(x, y, z)
    end
    @inbounds tensor[i, j, k] = sum
end

@kernel function dipolar_tensors_kernel_xy!(tensor, dx::Float64, dy::Float64, dz::Float64,
                                            Nx::Int64, Ny::Int64, Nz::Int64)
    nx, ny, nz = size(tensor)
    i, j, k = @index(Global, NTuple)

    sum = 0.0
    for p in (-Nx):Nx, q in (-Ny):Ny, s in (-Nz):Nz
        x = (i - 1 + p * nx) * dx
        y = (j - 1 + q * ny) * dy
        z = (k - 1 + s * nz) * dz
        sum += dipolar_tensor_xy(x, y, z)
    end
    @inbounds tensor[i, j, k] = sum
end

@kernel function dipolar_tensors_kernel_xz!(tensor, dx::Float64, dy::Float64, dz::Float64,
                                            Nx::Int64, Ny::Int64, Nz::Int64)
    nx, ny, nz = size(tensor)
    i, j, k = @index(Global, NTuple)

    sum = 0.0
    for p in (-Nx):Nx, q in (-Ny):Ny, s in (-Nz):Nz
        x = (i - 1 + p * nx) * dx
        y = (j - 1 + q * ny) * dy
        z = (k - 1 + s * nz) * dz
        sum += dipolar_tensor_xz(x, y, z)
    end
    @inbounds tensor[i, j, k] = sum
end

@kernel function dipolar_tensors_kernel_yz!(tensor, dx::Float64, dy::Float64, dz::Float64,
                                            Nx::Int64, Ny::Int64, Nz::Int64)
    nx, ny, nz = size(tensor)
    i, j, k = @index(Global, NTuple)

    sum = 0.0
    for p in (-Nx):Nx, q in (-Ny):Ny, s in (-Nz):Nz
        x = (i - 1 + p * nx) * dx
        y = (j - 1 + q * ny) * dy
        z = (k - 1 + s * nz) * dz
        sum += dipolar_tensor_yz(x, y, z)
    end
    @inbounds tensor[i, j, k] = sum
end

function compute_dipolar_tensors(tensor, kernel_fun, Nx, Ny, Nz, dx, dy, dz)
    kernel! = kernel_fun(get_backend(tensor), groupsize[])
    kernel!(tensor, dx, dy, dz, Nx, Ny, Nz; ndrange=size(tensor))
    KernelAbstractions.synchronize(get_backend(tensor))
    return nothing
end

@kernel function distribute_m_atomistic_kernel!(@Const(m), mx_pad, my_pad, mz_pad,
                                                @Const(mu_s))
    i, j, k = @index(Global, NTuple)
    I = @index(Global)

    p = 3 * I - 2
    @inbounds mx_pad[i, j, k] = m[p] * mu_s[I] * 1e20  # 1e20 = 4*pi*mu_0*(1e9)^3
    @inbounds my_pad[i, j, k] = m[p + 1] * mu_s[I] * 1e20
    @inbounds mz_pad[i, j, k] = m[p + 2] * mu_s[I] * 1e20
end

function distribute_m_atomistic(m, mx_pad, my_pad, mz_pad, mu_s::AbstractArray{T,1},
                                nx::Int64, ny::Int64, nz::Int64) where {T<:AbstractFloat}
    kernel! = distribute_m_atomistic_kernel!(default_backend[], groupsize[])
    kernel!(m, mx_pad, my_pad, mz_pad, mu_s; ndrange=(nx, ny, nz))
    KernelAbstractions.synchronize(default_backend[])
    return nothing
end

@kernel function collect_h_atomistic_kernel!(h, energy, @Const(m), @Const(hx), @Const(hy),
                                             @Const(hz), @Const(mu_s))
    i, j, k = @index(Global, NTuple)
    I = @index(Global)

    p = 3 * I - 2
    @inbounds h[p] = -hx[i, j, k]
    @inbounds h[p + 1] = -hy[i, j, k]
    @inbounds h[p + 2] = -hz[i, j, k]

    @inbounds mh = m[p] * h[p] + m[p + 1] * h[p + 1] + m[p + 2] * h[p + 2]
    @inbounds energy[I] = -0.5 * mu_s[I] * mh
end

function collect_h_atomistic_energy(h, energy, m, hx, hy, hz, mu_s::AbstractArray{T,1},
                                    nx::Int64, ny::Int64,
                                    nz::Int64) where {T<:AbstractFloat}
    kernel! = collect_h_atomistic_kernel!(default_backend[], groupsize[])
    kernel!(h, energy, m, hx, hy, hz, mu_s; ndrange=(nx, ny, nz))
    KernelAbstractions.synchronize(default_backend[])
    return nothing
end
