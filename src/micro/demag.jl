using FFTW
using LinearAlgebra

mutable struct Demag{T<:AbstractFloat} <: MicroEnergy
    nx_fft::Int64
    ny_fft::Int64
    nz_fft::Int64
    tensor_xx::AbstractArray{T,3}
    tensor_yy::AbstractArray{T,3}
    tensor_zz::AbstractArray{T,3}
    tensor_xy::AbstractArray{T,3}
    tensor_xz::AbstractArray{T,3}
    tensor_yz::AbstractArray{T,3}
    mx::AbstractArray{T,3}  #input for FFT
    my::AbstractArray{T,3}
    mz::AbstractArray{T,3}
    Mx::AbstractArray{Complex{T},3} #output for FFT
    My::AbstractArray{Complex{T},3}
    Mz::AbstractArray{Complex{T},3}
    Hx::AbstractArray{Complex{T},3}
    Hy::AbstractArray{Complex{T},3}
    Hz::AbstractArray{Complex{T},3}
    m_plan::Any
    h_plan::Any
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    name::String
end

function init_demag(sim::MicroSim, Nx::Int, Ny::Int, Nz::Int)
    mesh = sim.mesh
    max_size = max(mesh.dx, mesh.dy, mesh.dz)
    dx = Float64(mesh.dx / max_size)
    dy = Float64(mesh.dy / max_size)
    dz = Float64(mesh.dz / max_size)

    nx = mesh.nx
    ny = mesh.ny
    nz = mesh.nz

    cn = 10
    nx_fft = mesh.nx > cn ? 2 * mesh.nx : 2 * mesh.nx - 1
    ny_fft = mesh.ny > cn ? 2 * mesh.ny : 2 * mesh.ny - 1
    nz_fft = mesh.nz > cn ? 2 * mesh.nz : 2 * mesh.nz - 1

    T = single_precision.x ? Float32 : Float64
    mx_pad = KernelAbstractions.zeros(backend[], T, nx_fft, ny_fft, nz_fft)
    my_pad = KernelAbstractions.zeros(backend[], T, nx_fft, ny_fft, nz_fft)
    mz_pad = KernelAbstractions.zeros(backend[], T, nx_fft, ny_fft, nz_fft)
    plan = plan_rfft(mx_pad)

    tensor = KernelAbstractions.zeros(backend[], T, nx, ny, nz)

    #Nxx
    compute_demag_tensors(tensor, tensors_kernel_xx!, Nx, Ny, Nz, dx, dy, dz)
    fill_tensors(mx_pad, tensor, false, false, false)
    tensor_xx = real(plan * mx_pad)

    #Nyy
    compute_demag_tensors(tensor, tensors_kernel_yy!, Nx, Ny, Nz, dx, dy, dz)
    fill_tensors(mx_pad, tensor, false, false, false)
    tensor_yy = real(plan * mx_pad)

    #Nzz
    compute_demag_tensors(tensor, tensors_kernel_zz!, Nx, Ny, Nz, dx, dy, dz)
    fill_tensors(mx_pad, tensor, false, false, false)
    tensor_zz = real(plan * mx_pad)

    #Nxy
    compute_demag_tensors(tensor, tensors_kernel_xy!, Nx, Ny, Nz, dx, dy, dz)
    fill_tensors(mx_pad, tensor, true, true, false)
    tensor_xy = real(plan * mx_pad)

    #Nxz
    compute_demag_tensors(tensor, tensors_kernel_xz!, Nx, Ny, Nz, dx, dy, dz)
    fill_tensors(mx_pad, tensor, true, false, true)
    tensor_xz = real(plan * mx_pad)

    #Nyz
    compute_demag_tensors(tensor, tensors_kernel_yz!, Nx, Ny, Nz, dx, dy, dz)
    fill_tensors(mx_pad, tensor, false, true, true)
    tensor_yz = real(plan * mx_pad)

    lenx = (nx_fft % 2 > 0) ? nx : nx + 1

    Mx = KernelAbstractions.zeros(backend[], Complex{T}, lenx, ny_fft, nz_fft)
    My = KernelAbstractions.zeros(backend[], Complex{T}, lenx, ny_fft, nz_fft)
    Mz = KernelAbstractions.zeros(backend[], Complex{T}, lenx, ny_fft, nz_fft)
    Hx = KernelAbstractions.zeros(backend[], Complex{T}, lenx, ny_fft, nz_fft)
    Hy = KernelAbstractions.zeros(backend[], Complex{T}, lenx, ny_fft, nz_fft)
    Hz = KernelAbstractions.zeros(backend[], Complex{T}, lenx, ny_fft, nz_fft)

    m_plan = plan_rfft(mx_pad)
    h_plan = plan_irfft(Hx, nx_fft)

    field = KernelAbstractions.zeros(backend[], T, 3 * sim.n_total)
    energy = KernelAbstractions.zeros(backend[], T, sim.n_total)

    demag = Demag(nx_fft, ny_fft, nz_fft, tensor_xx, tensor_yy, tensor_zz, tensor_xy,
                  tensor_xz, tensor_yz, mx_pad, my_pad, mz_pad, Mx, My, Mz, Hx, Hy, Hz,
                  m_plan, h_plan, field, energy, "Demag")
    return demag
end

function effective_field(demag::Demag, sim::MicroSim, spin::AbstractArray{T,1},
                         t::Float64) where {T<:AbstractFloat}
    mesh = sim.mesh
    nx, ny, nz = mesh.nx, sim.mesh.ny, sim.mesh.nz

    fill!(demag.mx, 0)
    fill!(demag.my, 0)
    fill!(demag.mz, 0)

    distribute_m(spin, demag.mx, demag.my, demag.mz, sim.Ms, nx, ny, nz)

    #synchronize()
    mul!(demag.Mx, demag.m_plan, demag.mx)
    mul!(demag.My, demag.m_plan, demag.my)
    mul!(demag.Mz, demag.m_plan, demag.mz)

    demag.Hx .= demag.tensor_xx .* demag.Mx .+ demag.tensor_xy .* demag.My .+
                demag.tensor_xz .* demag.Mz
    demag.Hy .= demag.tensor_xy .* demag.Mx .+ demag.tensor_yy .* demag.My .+
                demag.tensor_yz .* demag.Mz
    demag.Hz .= demag.tensor_xz .* demag.Mx .+ demag.tensor_yz .* demag.My .+
                demag.tensor_zz .* demag.Mz
    #synchronize()

    mul!(demag.mx, demag.h_plan, demag.Hx)
    mul!(demag.my, demag.h_plan, demag.Hy)
    mul!(demag.mz, demag.h_plan, demag.Hz)

    collect_h_energy(demag.field, demag.energy, spin, demag.mx, demag.my, demag.mz, sim.Ms,
                     T(mesh.volume), nx, ny, nz)

    return nothing
end

function newell_f(x::Float64, y::Float64, z::Float64)::Float64
    x2 = x * x
    y2 = y * y
    z2 = z * z
    R = sqrt(x2 + y2 + z2)
    if R == 0.0
        return 0.0
    end

    f = 1.0 / 6 * (2 * x2 - y2 - z2) * R
    if x2 > 0
        f -= x * y * z * atan(y * z / (x * R))
    end

    if x2 + z2 > 0
        f += 0.5 * y * (z2 - x2) * asinh(y / (sqrt(x2 + z2)))
    end

    if x2 + y2 > 0
        f += 0.5 * z * (y2 - x2) * asinh(z / (sqrt(x2 + y2)))
    end
    return f
end

function newell_g(x::Float64, y::Float64, z::Float64)::Float64
    x2 = x * x
    y2 = y * y
    z2 = z * z

    R = sqrt(x2 + y2 + z2)
    if R == 0.0
        return 0.0
    end

    g = -1.0 / 3 * x * y * R

    if z2 > 0
        g -= 1.0 / 6 * z2 * z * atan(x * y / (z * R))
    end
    if y2 > 0
        g -= 0.5 * y2 * z * atan(x * z / (y * R))
    end
    if x2 > 0
        g -= 0.5 * x2 * z * atan(y * z / (x * R))
    end

    if x2 + y2 > 0
        g += x * y * z * asinh(z / (sqrt(x2 + y2)))
    end

    if y2 + z2 > 0
        g += 1.0 / 6 * y * (3 * z2 - y2) * asinh(x / (sqrt(y2 + z2)))
    end

    if x2 + z2 > 0
        g += 1.0 / 6 * x * (3 * z2 - x2) * asinh(y / (sqrt(x2 + z2)))
    end

    return g
end

#Numerical Micromagnetics: Finite Difference Methods, Jacques E. Miltat1 and Michael J. Donahue. Page 14.

function demag_tensor_xx(x::Float64, y::Float64, z::Float64, dx::Float64, dy::Float64,
                         dz::Float64)
    tensor = 8.0 * newell_f(x, y, z)

    tensor -= 4.0 * newell_f(x + dx, y, z)
    tensor -= 4.0 * newell_f(x - dx, y, z)
    tensor -= 4.0 * newell_f(x, y - dy, z)
    tensor -= 4.0 * newell_f(x, y + dy, z)
    tensor -= 4.0 * newell_f(x, y, z - dz)
    tensor -= 4.0 * newell_f(x, y, z + dz)

    tensor += 2.0 * newell_f(x + dx, y + dy, z)
    tensor += 2.0 * newell_f(x + dx, y - dy, z)
    tensor += 2.0 * newell_f(x - dx, y - dy, z)
    tensor += 2.0 * newell_f(x - dx, y + dy, z)
    tensor += 2.0 * newell_f(x + dx, y, z + dz)
    tensor += 2.0 * newell_f(x + dx, y, z - dz)
    tensor += 2.0 * newell_f(x - dx, y, z + dz)
    tensor += 2.0 * newell_f(x - dx, y, z - dz)
    tensor += 2.0 * newell_f(x, y - dy, z - dz)
    tensor += 2.0 * newell_f(x, y - dy, z + dz)
    tensor += 2.0 * newell_f(x, y + dy, z + dz)
    tensor += 2.0 * newell_f(x, y + dy, z - dz)

    tensor -= newell_f(x + dx, y + dy, z + dz)
    tensor -= newell_f(x + dx, y + dy, z - dz)
    tensor -= newell_f(x + dx, y - dy, z + dz)
    tensor -= newell_f(x + dx, y - dy, z - dz)
    tensor -= newell_f(x - dx, y + dy, z + dz)
    tensor -= newell_f(x - dx, y + dy, z - dz)
    tensor -= newell_f(x - dx, y - dy, z + dz)
    tensor -= newell_f(x - dx, y - dy, z - dz)

    return tensor / (4.0 * pi * dx * dy * dz)
end

@inline function demag_tensor_yy(x::Float64, y::Float64, z::Float64, dx::Float64,
                                 dy::Float64, dz::Float64)
    return demag_tensor_xx(y, x, z, dy, dx, dz)
end

@inline function demag_tensor_zz(x::Float64, y::Float64, z::Float64, dx::Float64,
                                 dy::Float64, dz::Float64)
    return demag_tensor_xx(z, y, x, dz, dy, dx)
end

function demag_tensor_xy(x::Float64, y::Float64, z::Float64, dx::Float64, dy::Float64,
                         dz::Float64)
    tensor = 8.0 * newell_g(x, y, z)

    tensor -= 4.0 * newell_g(x + dx, y, z)
    tensor -= 4.0 * newell_g(x - dx, y, z)
    tensor -= 4.0 * newell_g(x, y - dy, z)
    tensor -= 4.0 * newell_g(x, y + dy, z)
    tensor -= 4.0 * newell_g(x, y, z - dz)
    tensor -= 4.0 * newell_g(x, y, z + dz)

    tensor += 2.0 * newell_g(x + dx, y + dy, z)
    tensor += 2.0 * newell_g(x + dx, y - dy, z)
    tensor += 2.0 * newell_g(x - dx, y - dy, z)
    tensor += 2.0 * newell_g(x - dx, y + dy, z)
    tensor += 2.0 * newell_g(x + dx, y, z + dz)
    tensor += 2.0 * newell_g(x + dx, y, z - dz)
    tensor += 2.0 * newell_g(x - dx, y, z + dz)
    tensor += 2.0 * newell_g(x - dx, y, z - dz)
    tensor += 2.0 * newell_g(x, y - dy, z - dz)
    tensor += 2.0 * newell_g(x, y - dy, z + dz)
    tensor += 2.0 * newell_g(x, y + dy, z + dz)
    tensor += 2.0 * newell_g(x, y + dy, z - dz)

    tensor -= newell_g(x + dx, y + dy, z + dz)
    tensor -= newell_g(x + dx, y + dy, z - dz)
    tensor -= newell_g(x + dx, y - dy, z + dz)
    tensor -= newell_g(x + dx, y - dy, z - dz)
    tensor -= newell_g(x - dx, y + dy, z + dz)
    tensor -= newell_g(x - dx, y + dy, z - dz)
    tensor -= newell_g(x - dx, y - dy, z + dz)
    tensor -= newell_g(x - dx, y - dy, z - dz)

    return tensor / (4.0 * pi * dx * dy * dz)
end

@inline function demag_tensor_xz(x::Float64, y::Float64, z::Float64, dx::Float64,
                                 dy::Float64, dz::Float64)
    return demag_tensor_xy(x, z, y, dx, dz, dy)
end

@inline function demag_tensor_yz(x::Float64, y::Float64, z::Float64, dx::Float64,
                                 dy::Float64, dz::Float64)
    return demag_tensor_xy(y, z, x, dy, dz, dx)
end

@kernel function tensors_kernel_xx!(tensor, dx::Float64, dy::Float64, dz::Float64,
                                    Nx::Int64, Ny::Int64, Nz::Int64)
    nx, ny, nz = size(tensor)
    i, j, k = @index(Global, NTuple)

    sum = 0.0
    for p in (-Nx):Nx, q in (-Ny):Ny, s in (-Nz):Nz
        x = (i - 1 + p * nx) * dx
        y = (j - 1 + q * ny) * dy
        z = (k - 1 + s * nz) * dz
        sum += demag_tensor_xx(x, y, z, dx, dy, dz)
    end
    @inbounds tensor[i, j, k] = sum
end

@kernel function tensors_kernel_yy!(tensor, dx::Float64, dy::Float64, dz::Float64,
                                    Nx::Int64, Ny::Int64, Nz::Int64)
    nx, ny, nz = size(tensor)
    i, j, k = @index(Global, NTuple)

    sum = 0.0
    for p in (-Nx):Nx, q in (-Ny):Ny, s in (-Nz):Nz
        x = (i - 1 + p * nx) * dx
        y = (j - 1 + q * ny) * dy
        z = (k - 1 + s * nz) * dz
        sum += demag_tensor_yy(x, y, z, dx, dy, dz)
    end
    @inbounds tensor[i, j, k] = sum
end

@kernel function tensors_kernel_zz!(tensor, dx::Float64, dy::Float64, dz::Float64,
                                    Nx::Int64, Ny::Int64, Nz::Int64)
    nx, ny, nz = size(tensor)
    i, j, k = @index(Global, NTuple)

    sum = 0.0
    for p in (-Nx):Nx, q in (-Ny):Ny, s in (-Nz):Nz
        x = (i - 1 + p * nx) * dx
        y = (j - 1 + q * ny) * dy
        z = (k - 1 + s * nz) * dz
        sum += demag_tensor_zz(x, y, z, dx, dy, dz)
    end
    @inbounds tensor[i, j, k] = sum
end

@kernel function tensors_kernel_xy!(tensor, dx::Float64, dy::Float64, dz::Float64,
                                    Nx::Int64, Ny::Int64, Nz::Int64)
    nx, ny, nz = size(tensor)
    i, j, k = @index(Global, NTuple)

    sum = 0.0
    for p in (-Nx):Nx, q in (-Ny):Ny, s in (-Nz):Nz
        x = (i - 1 + p * nx) * dx
        y = (j - 1 + q * ny) * dy
        z = (k - 1 + s * nz) * dz
        sum += demag_tensor_xy(x, y, z, dx, dy, dz)
    end
    @inbounds tensor[i, j, k] = sum
end

@kernel function tensors_kernel_xz!(tensor, dx::Float64, dy::Float64, dz::Float64,
                                    Nx::Int64, Ny::Int64, Nz::Int64)
    nx, ny, nz = size(tensor)
    i, j, k = @index(Global, NTuple)

    sum = 0.0
    for p in (-Nx):Nx, q in (-Ny):Ny, s in (-Nz):Nz
        x = (i - 1 + p * nx) * dx
        y = (j - 1 + q * ny) * dy
        z = (k - 1 + s * nz) * dz
        sum += demag_tensor_xz(x, y, z, dx, dy, dz)
    end
    @inbounds tensor[i, j, k] = sum
end

@kernel function tensors_kernel_yz!(tensor, dx::Float64, dy::Float64, dz::Float64,
                                    Nx::Int64, Ny::Int64, Nz::Int64)
    nx, ny, nz = size(tensor)
    i, j, k = @index(Global, NTuple)

    sum = 0.0
    for p in (-Nx):Nx, q in (-Ny):Ny, s in (-Nz):Nz
        x = (i - 1 + p * nx) * dx
        y = (j - 1 + q * ny) * dy
        z = (k - 1 + s * nz) * dz
        sum += demag_tensor_yz(x, y, z, dx, dy, dz)
    end
    @inbounds tensor[i, j, k] = sum
end

function compute_demag_tensors(tensor, kernel_fun, Nx, Ny, Nz, dx, dy, dz)
    groupsize = 512
    kernel! = kernel_fun(get_backend(tensor), groupsize)
    kernel!(tensor, dx, dy, dz, Nx, Ny, Nz; ndrange=size(tensor))
    KernelAbstractions.synchronize(get_backend(tensor))
    return nothing
end

@kernel function fill_tensors_kernel!(long_tensor, tensor, tx::Bool, ty::Bool, tz::Bool)
    lnx, lny, lnz = size(long_tensor)
    nx, ny, nz = size(tensor)
    i, j, k = @index(Global, NTuple)

    if (lnx % 2 == 0 && i == nx + 1) ||
       (lny % 2 == 0 && j == ny + 1) ||
       (lnz % 2 == 0 && k == nz + 1)
    else
        x = (i <= nx) ? i : lnx - i + 2
        y = (j <= ny) ? j : lny - j + 2
        z = (k <= nz) ? k : lnz - k + 2
        sx = tx && (i > nx) ? -1 : 1
        sy = ty && (j > ny) ? -1 : 1
        sz = tz && (k > nz) ? -1 : 1
        @inbounds long_tensor[i, j, k] = sx * sy * sz * tensor[x, y, z]
    end
end

function fill_tensors(long_tensor, tensor, tx::Bool, ty::Bool, tz::Bool)
    groupsize = 512
    kernel! = fill_tensors_kernel!(get_backend(tensor), groupsize)
    kernel!(long_tensor, tensor, tx, ty, tz; ndrange=size(long_tensor))
    KernelAbstractions.synchronize(get_backend(tensor))
    return nothing
end

@kernel function distribute_m_kernel!(@Const(m), mx_pad, my_pad, mz_pad, @Const(Ms))
    i, j, k = @index(Global, NTuple)
    I = @index(Global)

    p = 3 * I - 2
    @inbounds mx_pad[i, j, k] = m[p] * Ms[I]
    @inbounds my_pad[i, j, k] = m[p + 1] * Ms[I]
    @inbounds mz_pad[i, j, k] = m[p + 2] * Ms[I]
end

function distribute_m(m, mx_pad, my_pad, mz_pad, Ms, nx::Int64, ny::Int64, nz::Int64)
    groupsize = 512
    kernel! = distribute_m_kernel!(backend[], groupsize)
    kernel!(m, mx_pad, my_pad, mz_pad, Ms; ndrange=(nx, ny, nz))
    KernelAbstractions.synchronize(backend[])
    return nothing
end

@kernel function collect_h_kernel!(h, energy, @Const(m), @Const(hx), @Const(hy), @Const(hz),
                                   @Const(Ms), volume::T) where {T<:AbstractFloat}
    i, j, k = @index(Global, NTuple)
    I = @index(Global)

    p = 3 * I - 2
    @inbounds h[p] = -1.0 * hx[i, j, k]
    @inbounds h[p + 1] = -1.0 * hy[i, j, k]
    @inbounds h[p + 2] = -1.0 * hz[i, j, k]

    @inbounds mh = m[p] * h[p] + m[p + 1] * h[p + 1] + m[p + 2] * h[p + 2]

    @inbounds energy[I] = -0.5 * mu_0 * Ms[I] * volume * mh
end

function collect_h_energy(h, energy, m, hx, hy, hz, Ms, volume::T, nx::Int64, ny::Int64,
                          nz::Int64) where {T<:AbstractFloat}
    groupsize = 512
    kernel! = collect_h_kernel!(backend[], groupsize)
    kernel!(h, energy, m, hx, hy, hz, Ms, volume; ndrange=(nx, ny, nz))
    KernelAbstractions.synchronize(backend[])
    return nothing
end
