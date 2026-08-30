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

    lenx = (nx_fft % 2 > 0) ? nx : nx + 1

    T = Float[]
    #one batched buffer for the 3 magnetization components (see micro/demag.jl)
    m_pad = create_zeros(nx_fft, ny_fft, nz_fft, 3)
    M_pad = create_zeros(Complex{T}, lenx, ny_fft, nz_fft, 3)
    if inplace_inverse(M_pad)
        #in-place c2r: the real field overwrites the spectrum buffer's own
        #memory, so no separate h_pad buffer is needed; 1/N is folded into the
        #packed tensors (see the comments in micro/demag.jl)
        tune_fftw_threads(m_pad, M_pad, nx_fft, make_inplace_plan)
        h_pad = reinterpret(T, M_pad)
        h_plan = make_inplace_plan(M_pad, nx_fft)
        tscale = T(inv(nx_fft * ny_fft * nz_fft))
    else
        h_pad = create_zeros(nx_fft, ny_fft, nz_fft, 3)
        h_plan = plan_irfft(M_pad, nx_fft, 1:3)
        tscale = one(T)
    end

    m_plan = plan_rfft(m_pad, 1:3)

    tensor = create_zeros(nx, ny, nz)
    mx_pad = create_zeros(nx_fft, ny_fft, nz_fft)
    plan = plan_rfft(mx_pad)

    #Nxx
    compute_dipolar_tensors(tensor, dipolar_tensors_kernel_xx!, Nx, Ny, Nz, dx, dy, dz)
    fill_tensors(mx_pad, tensor, false, false, false)
    tensor_xx = pack_demag_tensor(real(plan * mx_pad), 1, 1; tscale = tscale)

    #Nyy
    compute_dipolar_tensors(tensor, dipolar_tensors_kernel_yy!, Nx, Ny, Nz, dx, dy, dz)
    fill_tensors(mx_pad, tensor, false, false, false)
    tensor_yy = pack_demag_tensor(real(plan * mx_pad), 1, 1; tscale = tscale)

    #Nzz
    compute_dipolar_tensors(tensor, dipolar_tensors_kernel_zz!, Nx, Ny, Nz, dx, dy, dz)
    fill_tensors(mx_pad, tensor, false, false, false)
    tensor_zz = pack_demag_tensor(real(plan * mx_pad), 1, 1; tscale = tscale)

    #Nxy
    compute_dipolar_tensors(tensor, dipolar_tensors_kernel_xy!, Nx, Ny, Nz, dx, dy, dz)
    fill_tensors(mx_pad, tensor, true, true, false)
    tensor_xy = pack_demag_tensor(real(plan * mx_pad), -1, 1; tscale = tscale)

    #Nxz
    compute_dipolar_tensors(tensor, dipolar_tensors_kernel_xz!, Nx, Ny, Nz, dx, dy, dz)
    fill_tensors(mx_pad, tensor, true, false, true)
    tensor_xz = pack_demag_tensor(real(plan * mx_pad), 1, -1; tscale = tscale)

    #Nyz
    compute_dipolar_tensors(tensor, dipolar_tensors_kernel_yz!, Nx, Ny, Nz, dx, dy, dz)
    fill_tensors(mx_pad, tensor, false, true, true)
    tensor_yz = pack_demag_tensor(real(plan * mx_pad), -1, -1; tscale = tscale)

    field = create_zeros(3 * sim.n_total)
    energy = create_zeros(sim.n_total)

    demag = Demag(nx_fft, ny_fft, nz_fft, tensor_xx, tensor_yy, tensor_zz, tensor_xy,
                  tensor_xz, tensor_yz, m_pad, M_pad, h_pad, m_plan, h_plan,
                  field, energy, "Demag")
    return demag
end

function effective_field(demag::Demag, sim::AtomisticSim, spin::AbstractArray{T,1},
                         t::Float64) where {T<:AbstractFloat}
    mesh = sim.mesh
    nx, ny, nz = mesh.nx, sim.mesh.ny, sim.mesh.nz

    #see the comment in micro/demag.jl: the zero padding of m_pad stays zero
    distribute_m_atomistic(spin, demag.m_pad, sim.mu_s, nx, ny, nz)

    # synchronize() is not needed here: on the CPU backend KA kernels are
    # synchronous; on the CUDA backend both KA (distribute_m_atomistic) and
    # cuFFT (mul! with the rfft plan) execute on the default stream, so stream
    # ordering guarantees the kernel completes before the FFT begins.
    mul!(demag.M_pad, demag.m_plan, demag.m_pad)

    add_tensor_M(demag.M_pad, demag.tensor_xx, demag.tensor_yy, demag.tensor_zz,
                 demag.tensor_xy, demag.tensor_xz, demag.tensor_yz,
                 demag.ny_fft, demag.nz_fft)
    # synchronize() not needed: same default-stream ordering guarantee as above
    # applies between the KA kernel (add_tensor_M) and the cuFFT mul! calls.

    inv_transform!(demag.h_pad, demag.M_pad, demag.h_plan, demag.nx_fft)

    collect_h_atomistic_energy(demag.field, demag.energy, spin, demag.h_pad,
                               sim.mu_s, nx, ny, nz)

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
    return nothing
end

@kernel function distribute_m_atomistic_kernel!(@Const(m), m_pad, @Const(mu_s))
    i, j, k = @index(Global, NTuple)
    I = @index(Global)

    p = 3 * I - 2
    @inbounds m_pad[i, j, k, 1] = m[p] * mu_s[I] * 1e20  # 1e20 = 4*pi*mu_0*(1e9)^3
    @inbounds m_pad[i, j, k, 2] = m[p + 1] * mu_s[I] * 1e20
    @inbounds m_pad[i, j, k, 3] = m[p + 2] * mu_s[I] * 1e20
end

function distribute_m_atomistic(m, m_pad, mu_s::AbstractArray{T,1},
                                nx::Int64, ny::Int64, nz::Int64) where {T<:AbstractFloat}
    kernel! = distribute_m_atomistic_kernel!(default_backend[], groupsize[])
    kernel!(m, m_pad, mu_s; ndrange=(nx, ny, nz))
    return nothing
end

@kernel function collect_h_atomistic_kernel!(h, energy, @Const(m), @Const(h_pad),
                                             @Const(mu_s))
    i, j, k = @index(Global, NTuple)
    I = @index(Global)

    p = 3 * I - 2
    @inbounds h[p] = -h_pad[i, j, k, 1]
    @inbounds h[p + 1] = -h_pad[i, j, k, 2]
    @inbounds h[p + 2] = -h_pad[i, j, k, 3]

    @inbounds mh = m[p] * h[p] + m[p + 1] * h[p + 1] + m[p + 2] * h[p + 2]
    @inbounds energy[I] = -0.5 * mu_s[I] * mh
end

function collect_h_atomistic_energy(h, energy, m, h_pad, mu_s::AbstractArray{T,1},
                                    nx::Int64, ny::Int64,
                                    nz::Int64) where {T<:AbstractFloat}
    kernel! = collect_h_atomistic_kernel!(default_backend[], groupsize[])
    kernel!(h, energy, m, h_pad, mu_s; ndrange=(nx, ny, nz))
    return nothing
end
