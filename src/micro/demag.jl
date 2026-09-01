using FFTW
using LinearAlgebra

mutable struct Demag{T<:AbstractFloat} <: MicroEnergy
    nx_fft::Int64
    ny_fft::Int64
    nz_fft::Int64
    #Fourier-space demag tensors packed into the fundamental domain of the
    #(y,z) mirror parity group: (lenx, ny_fft÷2+1, nz_fft÷2+1), see
    #pack_demag_tensor below. Roughly 4x smaller than the full half-spectrum.
    tensor_xx::AbstractArray{T,3}
    tensor_yy::AbstractArray{T,3}
    tensor_zz::AbstractArray{T,3}
    tensor_xy::AbstractArray{T,3}
    tensor_xz::AbstractArray{T,3}
    tensor_yz::AbstractArray{T,3}
    m_pad::AbstractArray{T,4}           #padded m (nx_fft, ny_fft, nz_fft, 3), zero padding is kept zero
    M_pad::AbstractArray{Complex{T},4}  #spectrum, batched over the 3 components
    h_pad::AbstractArray{T,4}           #padded h, only [1:nx, 1:ny, 1:nz] is meaningful;
                                        #a reinterpret view of M_pad on backends with an
                                        #in-place c2r (CPU/CUDA), an own buffer otherwise
    m_plan::Any
    h_plan::Any
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    name::String
end

#The padded demag kernels constructed by fill_tensors have exact mirror
#parities on the cyclic (nx_fft, ny_fft, nz_fft) grid,
#    F(a, ny_fft-b, c) = ey * F(a,b,c),   F(a, b, nz_fft-c) = ez * F(a,b,c),
#with (ey, ez) = (+1,+1) for the diagonal components xx/yy/zz, (-1,+1) for xy,
#(+1,-1) for xz and (-1,-1) for yz (F is real because the kernels are even
#under the full point reflection).  The Fourier tensors are therefore stored
#only on the fundamental domain b < ny_fft÷2+1, c < nz_fft÷2+1 -- about 4x
#less memory than the full half-spectrum -- and add_tensor_M reconstructs the
#mirrored entries by index mapping.
#`full` (lenx, ny_fft, nz_fft) is the raw `real(plan * mx_pad)` output; it is
#used as scratch for the parity self-check and then discarded.
function pack_demag_tensor(full::AbstractArray{T,3}, ey::Int, ez::Int;
                           tscale::T=one(T)) where {T<:AbstractFloat}
    lenx, ny, nz = size(full)
    ny2 = ny ÷ 2 + 1
    nz2 = nz ÷ 2 + 1
    packed = similar(full, lenx, ny2, nz2)
    packed .= @view(full[:, 1:ny2, 1:nz2])

    #self-check: overwrite `full` with |full - parity-reconstructed value|
    scale = maximum(abs, full)
    kernel! = _tensor_parity_kernel!(get_backend(full), groupsize[])
    kernel!(full, packed, ey, ez, lenx, ny, nz, ny2, nz2; ndrange=length(full))
    vmax = maximum(full)
    #tolerance: FFT roundoff is ~1e-8 relative for Float64 and ~1e-5 (relative
    #to the dominant entries) for Float32; components that are physically
    #near-zero (e.g. yz of a monolayer, scale ~1e-17) are judged against an
    #absolute floor of the same O(1) kernel units instead
    tol = T == Float32 ? 1e-2 : 1e-9
    if vmax > tol * max(scale, one(T))
        @warn("demag tensor parity check failed: max |F - F_mirror| / max|F| = $(vmax / max(scale, one(T))); " *
              "the packed demag tensors would give wrong fields!")
    end
    tscale == one(T) || rmul!(packed, tscale)
    return packed
end

@kernel function _tensor_parity_kernel!(full, @Const(packed), ey::Int, ez::Int,
                                        lenx::Int, ny::Int, nz::Int, ny2::Int, nz2::Int)
    i = @index(Global)
    i0 = i - 1
    r, a = divrem(i0, lenx)
    c, b = divrem(r, ny)
    bb = b < ny2 ? b : ny - b
    cc = c < nz2 ? c : nz - c
    sy = b < ny2 ? 1 : ey
    sz = c < nz2 ? 1 : ez
    @inbounds expected = sy * sz * packed[a + 1 + lenx * (bb + ny2 * cc)]
    @inbounds full[i] = abs(full[i] - expected)
end

#Backends with a raw in-place c2r transform (CPU/FFTW and CUDA/cuFFT) write the
#real field directly into the spectrum buffer's own memory: the padded real
#output (2*lenx = nx_fft+2 real rows per component) exactly fits the complex
#input (lenx rows), so the separate h_pad allocation disappears.  Other
#backends fall back to an independent h_pad plus the public plan API.  The raw
#c2r transform is UNNORMALIZED, hence 1/N is folded into the packed tensors.
inplace_inverse(::AbstractArray) = false
inplace_inverse(A::Array) = true
inplace_inverse(A::Array{Complex{T}}) where {T} = T === Float64 || T === Float32

#raw FFTW in-place c2r plan, batched over the 3 components (howmany = 3)
mutable struct FFTWInplacePlan
    ptr::Ptr{Cvoid}
    destroy::Function
end

function make_inplace_plan(M_pad::Array{Complex{Float64},4}, nx_fft::Int)
    lenx, ny, nz = size(M_pad)[1:3]
    n = Cint[nz, ny, nx_fft]            #C order: the halved (first) dim goes last
    inembed = Cint[nz, ny, lenx]        #complex input dims
    onembed = Cint[nz, ny, 2 * lenx]    #padded real output dims
    ptr = ccall((:fftw_plan_many_dft_c2r, FFTW.libfftw3), Ptr{Cvoid},
                (Cint, Ptr{Cint}, Cint, Ptr{ComplexF64}, Ptr{Cint}, Cint, Cint,
                 Ptr{Float64}, Ptr{Cint}, Cint, Cint, Cuint),
                Cint(3), n, Cint(3), pointer(M_pad), inembed, Cint(1),
                Cint(lenx * ny * nz), pointer(M_pad), onembed, Cint(1),
                Cint(2 * lenx * ny * nz), Cuint(FFTW.MEASURE))
    ptr == C_NULL && error("fftw_plan_many_dft_c2r failed")
    p = FFTWInplacePlan(ptr, () -> ccall((:fftw_destroy_plan, FFTW.libfftw3), Cvoid,
                                         (Ptr{Cvoid},), ptr))
    finalizer(p) do pl
        pl.destroy()
    end
    return p
end

function make_inplace_plan(M_pad::Array{Complex{Float32},4}, nx_fft::Int)
    lenx, ny, nz = size(M_pad)[1:3]
    n = Cint[nz, ny, nx_fft]
    inembed = Cint[nz, ny, lenx]
    onembed = Cint[nz, ny, 2 * lenx]
    ptr = ccall((:fftwf_plan_many_dft_c2r, FFTW.libfftw3f), Ptr{Cvoid},
                (Cint, Ptr{Cint}, Cint, Ptr{ComplexF32}, Ptr{Cint}, Cint, Cint,
                 Ptr{Float32}, Ptr{Cint}, Cint, Cint, Cuint),
                Cint(3), n, Cint(3), pointer(M_pad), inembed, Cint(1),
                Cint(lenx * ny * nz), pointer(M_pad), onembed, Cint(1),
                Cint(2 * lenx * ny * nz), Cuint(FFTW.MEASURE))
    ptr == C_NULL && error("fftwf_plan_many_dft_c2r failed")
    p = FFTWInplacePlan(ptr, () -> ccall((:fftwf_destroy_plan, FFTW.libfftw3f), Cvoid,
                                         (Ptr{Cvoid},), ptr))
    finalizer(p) do pl
        pl.destroy()
    end
    return p
end

#execute the inverse: c2r with in == out (destroys M_pad, refilled every step)
function inv_transform!(h_pad, M_pad::Array{Complex{Float64},4}, h_plan::FFTWInplacePlan,
                        nx_fft::Int)
    ccall((:fftw_execute_dft_c2r, FFTW.libfftw3), Cvoid,
          (Ptr{Cvoid}, Ptr{Float64}, Ptr{Float64}),
          h_plan.ptr, pointer(M_pad), pointer(M_pad))
    return nothing
end

function inv_transform!(h_pad, M_pad::Array{Complex{Float32},4}, h_plan::FFTWInplacePlan,
                        nx_fft::Int)
    ccall((:fftwf_execute_dft_c2r, FFTW.libfftw3f), Cvoid,
          (Ptr{Cvoid}, Ptr{Float32}, Ptr{Float32}),
          h_plan.ptr, pointer(M_pad), pointer(M_pad))
    return nothing
end

#generic fallback for backends without a raw in-place c2r
inv_transform!(h_pad, M_pad, h_plan, nx_fft::Int) = mul!(h_pad, h_plan, M_pad)

# Note: the size of the demag tensors is reduced to ~nx*ny*nz by packing them
# into the (y,z) parity fundamental domain, see pack_demag_tensor.
# FIXME: add the real pbc (current imeplenation is macro pbc)
#pbc1d_axis != 0 turns the macro image sum into a true 1D-periodic kernel: the
#nonzero one of Nx/Ny/Nz is the image count Ic along the periodic axis and the
#analytic far-field tail (demag_pbc1d.jl) is added to every quadrant tensor
function init_demag(sim::MicroSim, Nx::Int, Ny::Int, Nz::Int; pbc1d_axis::Int=0)
    mesh = sim.mesh
    max_size = max(mesh.dx, mesh.dy, mesh.dz)
    dx = Float64(mesh.dx / max_size)
    dy = Float64(mesh.dy / max_size)
    dz = Float64(mesh.dz / max_size)

    nx = mesh.nx
    ny = mesh.ny
    nz = mesh.nz
    pbc1d_n = pbc1d_axis == 1 ? Nx : pbc1d_axis == 2 ? Ny : Nz

    cn = 3
    nx_fft = mesh.nx > cn ? 2 * mesh.nx : 2 * mesh.nx - 1
    ny_fft = mesh.ny > cn ? 2 * mesh.ny : 2 * mesh.ny - 1
    nz_fft = mesh.nz > cn ? 2 * mesh.nz : 2 * mesh.nz - 1

    lenx = (nx_fft % 2 > 0) ? nx : nx + 1

    T = Float[]
    #one batched buffer for the 3 magnetization components (avoids 3 plan calls
    #per FFT stage and enables batched transforms on GPU backends)
    m_pad = create_zeros(nx_fft, ny_fft, nz_fft, 3)
    M_pad = create_zeros(Complex{T}, lenx, ny_fft, nz_fft, 3)
    if inplace_inverse(M_pad)
        #in-place c2r: the real field overwrites the spectrum buffer's own
        #memory (nx_fft+2 real rows fit into lenx complex rows per component),
        #so no separate h_pad buffer is needed; the raw transform is
        #unnormalized, hence 1/N is folded into the packed tensors below.
        #The thread count must be settled (tuned with throw-away plans) BEFORE
        #the final in-place plan is created.
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
    compute_demag_tensors(tensor, tensors_kernel_xx!, Nx, Ny, Nz, dx, dy, dz)
    pbc1d_axis != 0 && pbc1d_add_tail!(tensor, pbc1d_axis, 1, pbc1d_n, dx, dy, dz, nx, ny, nz)
    pbc1d_axis != 0 && pbc1d_dc_fix!(tensor, pbc1d_axis, 1, dx, dy, dz, nx, ny, nz)
    fill_tensors(mx_pad, tensor, false, false, false)
    tensor_xx = pack_demag_tensor(real(plan * mx_pad), 1, 1; tscale = tscale)

    #Nyy
    compute_demag_tensors(tensor, tensors_kernel_yy!, Nx, Ny, Nz, dx, dy, dz)
    pbc1d_axis != 0 && pbc1d_add_tail!(tensor, pbc1d_axis, 2, pbc1d_n, dx, dy, dz, nx, ny, nz)
    pbc1d_axis != 0 && pbc1d_dc_fix!(tensor, pbc1d_axis, 2, dx, dy, dz, nx, ny, nz)
    fill_tensors(mx_pad, tensor, false, false, false)
    tensor_yy = pack_demag_tensor(real(plan * mx_pad), 1, 1; tscale = tscale)

    #Nzz
    compute_demag_tensors(tensor, tensors_kernel_zz!, Nx, Ny, Nz, dx, dy, dz)
    pbc1d_axis != 0 && pbc1d_add_tail!(tensor, pbc1d_axis, 3, pbc1d_n, dx, dy, dz, nx, ny, nz)
    pbc1d_axis != 0 && pbc1d_dc_fix!(tensor, pbc1d_axis, 3, dx, dy, dz, nx, ny, nz)
    fill_tensors(mx_pad, tensor, false, false, false)
    tensor_zz = pack_demag_tensor(real(plan * mx_pad), 1, 1; tscale = tscale)

    #Nxy
    compute_demag_tensors(tensor, tensors_kernel_xy!, Nx, Ny, Nz, dx, dy, dz)
    pbc1d_axis != 0 && pbc1d_add_tail!(tensor, pbc1d_axis, 4, pbc1d_n, dx, dy, dz, nx, ny, nz)
    pbc1d_axis != 0 && pbc1d_dc_fix!(tensor, pbc1d_axis, 4, dx, dy, dz, nx, ny, nz)
    fill_tensors(mx_pad, tensor, true, true, false)
    tensor_xy = pack_demag_tensor(real(plan * mx_pad), -1, 1; tscale = tscale)

    #Nxz
    compute_demag_tensors(tensor, tensors_kernel_xz!, Nx, Ny, Nz, dx, dy, dz)
    pbc1d_axis != 0 && pbc1d_add_tail!(tensor, pbc1d_axis, 5, pbc1d_n, dx, dy, dz, nx, ny, nz)
    pbc1d_axis != 0 && pbc1d_dc_fix!(tensor, pbc1d_axis, 5, dx, dy, dz, nx, ny, nz)
    fill_tensors(mx_pad, tensor, true, false, true)
    tensor_xz = pack_demag_tensor(real(plan * mx_pad), 1, -1; tscale = tscale)

    #Nyz
    compute_demag_tensors(tensor, tensors_kernel_yz!, Nx, Ny, Nz, dx, dy, dz)
    pbc1d_axis != 0 && pbc1d_add_tail!(tensor, pbc1d_axis, 6, pbc1d_n, dx, dy, dz, nx, ny, nz)
    pbc1d_axis != 0 && pbc1d_dc_fix!(tensor, pbc1d_axis, 6, dx, dy, dz, nx, ny, nz)
    fill_tensors(mx_pad, tensor, false, true, true)
    tensor_yz = pack_demag_tensor(real(plan * mx_pad), -1, -1; tscale = tscale)

    field = create_zeros(3 * sim.n_total)
    energy = create_zeros(sim.n_total)

    demag = Demag(nx_fft, ny_fft, nz_fft, tensor_xx, tensor_yy, tensor_zz, tensor_xy,
                  tensor_xz, tensor_yz, m_pad, M_pad, h_pad, m_plan, h_plan,
                  field, energy, "Demag")
    return demag
end

function effective_field(demag::Demag, sim::MicroSim, spin::AbstractArray{T,1}, t::Float64;
                         output=nothing) where {T<:AbstractFloat}
    mesh = sim.mesh
    nx, ny, nz = mesh.nx, sim.mesh.ny, sim.mesh.nz

    # The zero padding of m_pad is zeroed once at init and never written again:
    # the forward r2c mul! preserves its input and distribute_m fully overwrites
    # the mesh region at every call, while the padded tail stays zero.
    distribute_m(spin, demag.m_pad, sim.mu0_Ms, nx, ny, nz)

    # synchronize() is not needed here: on the CPU backend KA kernels are
    # synchronous; on the CUDA backend both KA (distribute_m) and cuFFT (mul!
    # with the rfft plan) execute on the default stream, so stream ordering
    # guarantees the kernel completes before the FFT begins.
    mul!(demag.M_pad, demag.m_plan, demag.m_pad)

    add_tensor_M(demag.M_pad, demag.tensor_xx, demag.tensor_yy, demag.tensor_zz,
                 demag.tensor_xy, demag.tensor_xz, demag.tensor_yz,
                 demag.ny_fft, demag.nz_fft)

    #in-place c2r (CPU/CUDA): the real field overwrites M_pad's memory and h_pad
    #is a reinterpret view of it; other backends: out-of-place into own buffer
    inv_transform!(demag.h_pad, demag.M_pad, demag.h_plan, demag.nx_fft)

    heff = output == nothing ? demag.field : output

    collect_h_energy(heff, demag.energy, spin, demag.h_pad, sim.mu0_Ms,
                     T(mesh.volume), nx, ny, nz)

    return nothing
end

#FFTW only parallelizes plans that are created after set_num_threads, and the
#optimal thread count depends strongly on the grid shape (e.g. very thin grids
#can be several times SLOWER with many threads), so we time the REAL pipeline
#(forward mul! + in-place c2r via `make_ip_plan`) for both candidates and keep
#the faster one. The timing plans are thrown away; the caller creates the final
#plans afterwards with the winning thread count. Non-CPU backends (cuFFT etc.)
#ignore FFTW and are left untouched.
function tune_fftw_threads(m_pad::AbstractArray{T,4}, M_pad::AbstractArray{Complex{T},4},
                           nx_fft::Int, make_ip_plan::Function) where {T<:AbstractFloat}
    default_backend[] != CPU() && return nothing
    if length(m_pad) < 2^16 #threading cannot pay off on small problems
        FFTW.set_num_threads(1)
        return nothing
    end

    maxt = max(Sys.CPU_THREADS, 1)
    best, best_time = 1, Inf
    for nt in unique((1, maxt))
        FFTW.set_num_threads(nt)
        p = plan_rfft(m_pad, 1:3)
        ip = make_ip_plan(M_pad, nx_fft)
        #several rounds, keep the fastest: a single noisy sample on a busy
        #machine can otherwise select a catastrophically bad thread count
        t = Inf
        for _ in 1:4
            mul!(M_pad, p, m_pad)        #warmup + refill (r2c preserves m_pad)
            inv_transform!(M_pad, M_pad, ip, nx_fft)  #warmup (destroys M_pad)
            t = min(t, @elapsed begin
                mul!(M_pad, p, m_pad)
                inv_transform!(M_pad, M_pad, ip, nx_fft)
            end)
        end
        finalize(ip)
        if t < best_time
            best, best_time = nt, t
        end
    end
    FFTW.set_num_threads(best)
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
    R = sqrt(x * x + y * y + z * z) / max(dx, dy, dz)
    if R > 60 #use the dipolar approximation
        return dipolar_tensor_xx(x, y, z) * dx * dy * dz / (4 * pi)
    end

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
    R = sqrt(x * x + y * y + z * z) / max(dx, dy, dz)
    if R > 60 #use the dipolar approximation
        return dipolar_tensor_xy(x, y, z) * dx * dy * dz / (4 * pi)
    end

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
    kernel! = kernel_fun(get_backend(tensor), groupsize[])
    kernel!(tensor, dx, dy, dz, Nx, Ny, Nz; ndrange=size(tensor))
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
    kernel! = fill_tensors_kernel!(get_backend(tensor))
    kernel!(long_tensor, tensor, tx, ty, tz; ndrange=size(long_tensor))
    return nothing
end

@kernel function distribute_m_kernel!(@Const(m), m_pad, Ms)
    i, j, k = @index(Global, NTuple)
    I = @index(Global)

    p = 3 * I - 2
    @inbounds m_pad[i, j, k, 1] = m[p] * Ms[I]
    @inbounds m_pad[i, j, k, 2] = m[p + 1] * Ms[I]
    @inbounds m_pad[i, j, k, 3] = m[p + 2] * Ms[I]
end

function distribute_m(m, m_pad, Ms, nx::Int64, ny::Int64, nz::Int64)
    kernel! = distribute_m_kernel!(default_backend[])
    kernel!(m, m_pad, Ms; ndrange=(nx, ny, nz))
    return nothing
end

@kernel function collect_h_kernel!(h, energy, @Const(m), @Const(h_pad),
                                   mu0_Ms, volume::T) where {T<:AbstractFloat}
    i, j, k = @index(Global, NTuple)
    I = @index(Global)

    p = 3 * I - 2
    factor::T = -1.0 / mu_0
    @inbounds h[p] = factor * h_pad[i, j, k, 1]
    @inbounds h[p + 1] = factor * h_pad[i, j, k, 2]
    @inbounds h[p + 2] = factor * h_pad[i, j, k, 3]

    @inbounds mh = m[p] * h[p] + m[p + 1] * h[p + 1] + m[p + 2] * h[p + 2]

    @inbounds energy[I] = -0.5 * mu0_Ms[I] * volume * mh
end

function collect_h_energy(h, energy, m, h_pad, Ms, volume::T, nx::Int64, ny::Int64,
                          nz::Int64) where {T<:AbstractFloat}
    kernel! = collect_h_kernel!(default_backend[])
    kernel!(h, energy, m, h_pad, Ms, volume; ndrange=(nx, ny, nz))
    return nothing
end

# H, M and all tensors are 3D array, but we still can use 1d index here.
# Hx .= tensor_xx.*Mx .+ tensor_xy.*My .+  tensor_xz.*Mz
# Hy .= tensor_xy.*Mx .+ tensor_yy.*My .+  tensor_yz.*Mz
# Hz .= tensor_xz.*Mx .+ tensor_yz.*My .+  tensor_zz.*Mz
# we use M_pad to store Hx, Hy and Hz; the spectra of the three components are
# consecutive along the last (batch) dimension, so linear indices i, i + nspec
# and i + 2*nspec address the same frequency in Mx, My and Mz.
# The six tensors are stored packed on the (y,z) parity fundamental domain
# (see pack_demag_tensor): for a frequency (a, b, c) with b > ny2 or c > nz2
# the tensor entry is reconstructed from the root (a, ny+2-b, nz+2-c).  The
# diagonal components are even in both axes (no sign), xy is odd in y, xz odd
# in z and yz odd in both, giving the sign factors sy, sz below.  The kernel
# runs on a 3D ndrange (lenx, ny, nz) so that the frequency components are
# available directly (no integer division in the hot loop).
@kernel function add_tensor_M_kernel!(M_pad, @Const(tensor_xx), @Const(tensor_yy),
                                      @Const(tensor_zz), @Const(tensor_xy),
                                      @Const(tensor_xz), @Const(tensor_yz),
                                      lenx::Int64, ny::Int64, nz::Int64,
                                      ny2::Int64, nz2::Int64, nspec::Int64)
    ia, ib, ic = @index(Global, NTuple)
    bb = ib <= ny2 ? ib : ny + 2 - ib
    cc = ic <= nz2 ? ic : nz + 2 - ic
    sy = ib <= ny2 ? 1 : -1
    sz = ic <= nz2 ? 1 : -1
    i = ia + (ib - 1) * lenx + (ic - 1) * (lenx * ny)
    p = ia + (bb - 1) * lenx + (cc - 1) * (lenx * ny2)

    @inbounds Mx = M_pad[i]
    @inbounds My = M_pad[i + nspec]
    @inbounds Mz = M_pad[i + 2 * nspec]

    @inbounds xx = tensor_xx[p]
    @inbounds yy = tensor_yy[p]
    @inbounds zz = tensor_zz[p]
    @inbounds xy = sy * tensor_xy[p]
    @inbounds xz = sz * tensor_xz[p]
    @inbounds yz = sy * sz * tensor_yz[p]

    @inbounds Hx = xx * Mx + xy * My + xz * Mz
    @inbounds Hy = xy * Mx + yy * My + yz * Mz
    @inbounds Hz = xz * Mx + yz * My + zz * Mz

    @inbounds M_pad[i] = Hx
    @inbounds M_pad[i + nspec] = Hy
    @inbounds M_pad[i + 2 * nspec] = Hz
end

function add_tensor_M(M_pad, tensor_xx, tensor_yy, tensor_zz, tensor_xy, tensor_xz,
                      tensor_yz, ny_fft::Int64, nz_fft::Int64)
    kernel! = add_tensor_M_kernel!(default_backend[], groupsize[])
    lenx, ny2, nz2 = size(tensor_xx)
    nspec = lenx * ny_fft * nz_fft
    kernel!(M_pad, tensor_xx, tensor_yy, tensor_zz, tensor_xy, tensor_xz, tensor_yz,
            lenx, ny_fft, nz_fft, ny2, nz2, nspec; ndrange=(lenx, ny_fft, nz_fft))
    return nothing
end
