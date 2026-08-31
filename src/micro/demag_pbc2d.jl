#True 2D-periodic demag (x, y periodic; z open) for the FD grid.
#
#Method (Wang et al., Comp. Mater. Sci. 49 (2010) 84): the in-plane-
#periodized Newell kernel with the explicit images |p|<=Ic, |q|<=Jc plus the
#analytic far-field tail N_inf (the continuum integral of the dipolar kernel
#over the region outside the critical rectangle, closed forms Eq. 13; residual
#error ~ gamma^-3/2, estimable), transformed per mode by the in-plane DFT.
#The solve: 2D batched FFTs (no padding, no z-transform) + the O(nz^2) layer
#convolution + the 2D inverse FFT.
#The k=0 column (in-plane uniform) is analytic, not truncated: a uniform
#in-plane magnetization produces no field, and a uniform M_z is demagnetized
#by its own layer only (the +-M_z face-pair field is confined between its
#faces), i.e. H_z = -M_z per layer.
mutable struct DemagPBC2D{T<:AbstractFloat} <: MicroEnergy
    nx::Int64
    ny::Int64
    nz::Int64
    kernel::AbstractArray{Complex{T},4}   #(6, lenx, ny, 2nz-1) mode-space layer kernels
    m_pad::AbstractArray{T,4}             #(nx, ny, nz, 3) the periodic cell itself
    M_pad::AbstractArray{Complex{T},4}    #(lenx, ny, nz, 3) 2D spectrum
    M_out::AbstractArray{Complex{T},4}    #the layer-convolution output (separate
                                          #buffer: the convolution reads all layers
                                          #of M_pad, so writing in place would race)
    h_pad::AbstractArray{T,4}             #(nx, ny, nz, 3)
    m_plan::Any
    h_plan::Any
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    name::String
end

#F functions: the closed-form quarter-plane integrals of the dipolar kernel
#(the paper's Eq. 13); abc = the cell volume, Tx/Ty = the periods
function make_F2d(abc::T, Tx::T, Ty::T) where {T<:AbstractFloat}
    Fxx(p, q, z) = abc / (4pi * Tx * Ty) * p * q / ((p^2 + z^2) * sqrt(p^2 + q^2 + z^2))
    Fyy(p, q, z) = abc / (4pi * Tx * Ty) * p * q / ((q^2 + z^2) * sqrt(p^2 + q^2 + z^2))
    Fzz(p, q, z) = -abc / (4pi * Tx * Ty) * p * q * (p^2 + q^2 + 2z^2) /
        ((p^2 + z^2) * (q^2 + z^2) * sqrt(p^2 + q^2 + z^2))
    Fxy(p, q, z) = -abc / (4pi * Tx * Ty) / sqrt(p^2 + q^2 + z^2)
    Fxz(p, q, z) = abc / (4pi * Tx * Ty) * q * z / ((p^2 + z^2) * sqrt(p^2 + q^2 + z^2))
    Fyz(p, q, z) = abc / (4pi * Tx * Ty) * p * z / ((q^2 + z^2) * sqrt(p^2 + q^2 + z^2))
    return Fxx, Fyy, Fzz, Fxy, Fxz, Fyz
end

#N_inf: the far-field tail = the continuum integral of the dipolar kernel over
#the region outside the critical rectangle (the paper's Eq. 10-12); the four
#F terms combine the two semi-infinite x-strips, the two y-strips and the
#corners with the inclusion-exclusion signs
function Ninf2d(comp::Int, x::T, y::T, z::T, Xc::T, Yc::T, F) where {T<:AbstractFloat}
    Fxx, Fyy, Fzz, Fxy, Fxz, Fyz = F
    p1 = x + Xc; p2 = x - Xc
    q1 = y - Yc; q2 = y + Yc
    if comp == 1
        return Fxx(p1, q1, z) + Fxx(p2, q2, z) - Fxx(p1, q2, z) - Fxx(p2, q1, z)
    elseif comp == 2
        return Fyy(p1, q1, z) + Fyy(p2, q2, z) - Fyy(p1, q2, z) - Fyy(p2, q1, z)
    elseif comp == 3
        return Fzz(p1, q1, z) + Fzz(p2, q2, z) - Fzz(p1, q2, z) - Fzz(p2, q1, z)
    elseif comp == 4
        return Fxy(p1, q1, z) + Fxy(p2, q2, z) - Fxy(p1, q2, z) - Fxy(p2, q1, z)
    elseif comp == 5
        return Fxz(p1, q1, z) + Fxz(p2, q2, z) - Fxz(p1, q2, z) - Fxz(p2, q1, z)
    else
        return Fyz(p1, q1, z) + Fyz(p2, q2, z) - Fyz(p1, q2, z) - Fyz(p2, q1, z)
    end
end

function init_demag_pbc2d(sim::MicroSim; Ic::Int=2, Jc::Int=2)
    mesh = sim.mesh
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
    dx = Float64(mesh.dx); dy = Float64(mesh.dy); dz = Float64(mesh.dz)
    (nx < 2 || ny < 2) && error("pbc2d requires nx >= 2 and ny >= 2")
    T = Float[]
    lenx = nx ÷ 2 + 1
    m_pad = create_zeros(nx, ny, nz, 3)
    M_pad = create_zeros(Complex{T}, lenx, ny, nz, 3)
    M_out = create_zeros(Complex{T}, lenx, ny, nz, 3)
    h_pad = create_zeros(nx, ny, nz, 3)

    if default_backend[] == CPU()
        FFTW.set_num_threads(max(Sys.CPU_THREADS, 1))
    end
    m_plan = plan_rfft(m_pad, (1, 2))
    h_plan = plan_irfft(M_pad, nx, (1, 2))

    #---- kernel construction (CPU; one-time) ----
    abc = dx * dy * dz
    Tx = nx * dx; Ty = ny * dy
    Xc = (Ic + 0.5) * Tx; Yc = (Jc + 0.5) * Ty
    F = make_F2d(abc, Tx, Ty)
    Kreal = zeros(Float64, 6, nx, ny, 2nz - 1)   #comp, di, dj, n(+nz); built in F64
    for n in -(nz-1):(nz-1), dj in 0:(ny-1), di in 0:(nx-1)
        x = di * dx; y = dj * dy; z = n * dz
        #explicit images; the tensor components at the image displacement (the
        #axis permutation lives inside demag_tensor_ab, same as the macro path)
        for p in -Ic:Ic, q in -Jc:Jc
            xe = (di + p*nx) * dx; ye = (dj + q*ny) * dy
            Kreal[1, di+1, dj+1, n+nz] += demag_tensor_xx(xe, ye, z, dx, dy, dz)
            Kreal[2, di+1, dj+1, n+nz] += demag_tensor_xy(xe, ye, z, dx, dy, dz)
            Kreal[3, di+1, dj+1, n+nz] += demag_tensor_xz(xe, ye, z, dx, dy, dz)
            Kreal[4, di+1, dj+1, n+nz] += demag_tensor_yy(xe, ye, z, dx, dy, dz)
            Kreal[5, di+1, dj+1, n+nz] += demag_tensor_yz(xe, ye, z, dx, dy, dz)
            Kreal[6, di+1, dj+1, n+nz] += demag_tensor_zz(xe, ye, z, dx, dy, dz)
        end
        #analytic far-field tail
        Kreal[1, di+1, dj+1, n+nz] += Ninf2d(1, x, y, z, Xc, Yc, F)
        Kreal[2, di+1, dj+1, n+nz] += Ninf2d(2, x, y, z, Xc, Yc, F)
        Kreal[3, di+1, dj+1, n+nz] += Ninf2d(3, x, y, z, Xc, Yc, F)
        Kreal[4, di+1, dj+1, n+nz] += Ninf2d(4, x, y, z, Xc, Yc, F)
        Kreal[5, di+1, dj+1, n+nz] += Ninf2d(5, x, y, z, Xc, Yc, F)
        Kreal[6, di+1, dj+1, n+nz] += Ninf2d(6, x, y, z, Xc, Yc, F)
    end
    #in-plane DFT per (n, comp)
    kernel_f64 = zeros(ComplexF64, 6, lenx, ny, 2nz - 1)
    for nn in 1:(2nz-1), c in 1:6
        kernel_f64[c, :, :, nn] .= rfft(Kreal[c, :, :, nn], (1, 2))
    end
    #analytic k=0 column (the whole column): the in-plane-uniform modes see only
    #the per-layer self demag (the +-M_z face-pair field is confined between its
    #own faces; inter-layer coupling is exactly zero) -> H_z = -M_z, H_in = 0
    kernel_f64[:, 1, 1, :] .= 0.0
    kernel_f64[6, 1, 1, nz] = 1.0
    kernel = create_zeros(Complex{T}, 6, lenx, ny, 2nz - 1)
    copyto!(kernel, kernel_f64)

    field = create_zeros(3 * sim.n_total)
    energy = create_zeros(sim.n_total)
    return DemagPBC2D(nx, ny, nz, kernel, m_pad, M_pad, M_out, h_pad, m_plan, h_plan,
                      field, energy, "Demag")
end

#layer convolution per mode: H_i(f, m) = sum_{m'} N_ij(f, m-m') M_j(f, m');
#the tensor is symmetric, so the transposed products reuse the same kernel
#slots (N_yx = kernel[2] with Mx, N_zx = kernel[3] with Mx, N_zy = kernel[5] with My)
@kernel function pbc2d_conv_kernel!(M_out, @Const(M_pad), @Const(kernel), nz::Int64)
    a, b, m = @index(Global, NTuple)
    h1 = zero(M_pad[a, b, m, 1]); h2 = h1; h3 = h1
    for m2 in 1:nz
        n = (m - m2) + nz
        mx = M_pad[a, b, m2, 1]; my = M_pad[a, b, m2, 2]; mz = M_pad[a, b, m2, 3]
        k1 = kernel[1, a, b, n]; k2 = kernel[2, a, b, n]; k3 = kernel[3, a, b, n]
        k4 = kernel[4, a, b, n]; k5 = kernel[5, a, b, n]; k6 = kernel[6, a, b, n]
        h1 += k1 * mx + k2 * my + k3 * mz
        h2 += k2 * mx + k4 * my + k5 * mz
        h3 += k3 * mx + k5 * my + k6 * mz
    end
    @inbounds M_out[a, b, m, 1] = h1
    @inbounds M_out[a, b, m, 2] = h2
    @inbounds M_out[a, b, m, 3] = h3
end

function effective_field(demag::DemagPBC2D, sim::MicroSim, spin::AbstractArray{T,1},
                         t::Float64; output=nothing) where {T<:AbstractFloat}
    nx, ny, nz = demag.nx, demag.ny, demag.nz
    mesh = sim.mesh
    distribute_m(spin, demag.m_pad, sim.mu0_Ms, nx, ny, nz)
    mul!(demag.M_pad, demag.m_plan, demag.m_pad)
    kernel! = pbc2d_conv_kernel!(default_backend[], groupsize[])
    kernel!(demag.M_out, demag.M_pad, demag.kernel, nz;
            ndrange=(nx ÷ 2 + 1, ny, nz))
    mul!(demag.h_pad, demag.h_plan, demag.M_out)
    heff = output == nothing ? demag.field : output
    collect_h_energy(heff, demag.energy, spin, demag.h_pad, sim.mu0_Ms,
                     T(mesh.volume), nx, ny, nz)
    return nothing
end
