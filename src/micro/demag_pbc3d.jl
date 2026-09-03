#True 3D-periodic demag (tin-foil convention) for the FD grid.
#
#Instead of summing truncated dipolar images in real space (macro PBC, whose
#uniform mode is conditionally convergent and shape dependent), the periodic
#magnetostatic problem is solved directly in Fourier space: for every non-zero
#mode the kernel is the longitudinal projector
#    N_ij(k) = +k_i k_j / |k|^2
#(the sign convention of the open path: collect_h_energy applies -1/mu0),
#while the k = 0 (uniform) mode produces no field by convention (<H> = 0).
#The cell shape factor prod(sinc^2(pi*f/n)) makes this the exact periodic
#field of the BAND-LIMITED (trigonometric) reconstruction of the magnetization.
#The macro-PBC image sum converges instead to the top-hat (Newell) periodic
#kernel, whose longitudinal x-modes pick up the aliasing sum
#sum_m sinc^2(pi*(f/n + m)) = 1; the two discretizations share the same
#continuum limit and differ per mode by O(1 - prod(sinc^2)), shrinking ~dz^2
#under grid refinement.
#No zero padding is needed (circular = periodic convolution), so the FFTs are
#~8x smaller than the open-boundary path and no kernel tensors are stored.
mutable struct DemagPBC3D{T<:AbstractFloat} <: MicroEnergy
    nx::Int64
    ny::Int64
    nz::Int64
    m_pad::AbstractArray{T,4}           #(nx, ny, nz, 3), the periodic cell itself
    M_pad::AbstractArray{Complex{T},4}  #(nx÷2+1, ny, nz, 3) spectrum
    h_pad::AbstractArray{T,4}           #reinterpret view of M_pad (in-place c2r) or own buffer
    tscale::T                           #1/N for the raw unnormalized in-place c2r, 1 otherwise
    m_plan::Any
    h_plan::Any
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    name::String
end

function init_demag_pbc3d(sim::MicroSim)
    mesh = sim.mesh
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
    T = eltype(sim.spin)
    m_pad = create_zeros(nx, ny, nz, 3)
    M_pad = create_zeros(Complex{T}, nx ÷ 2 + 1, ny, nz, 3)
    if inplace_inverse(M_pad)
        tune_fftw_threads(m_pad, M_pad, nx, make_inplace_plan)
        h_pad = reinterpret(T, M_pad)
        h_plan = make_inplace_plan(M_pad, nx)
        tscale = T(inv(nx * ny * nz))
    else
        h_pad = create_zeros(nx, ny, nz, 3)
        h_plan = plan_irfft(M_pad, nx, 1:3)
        tscale = one(T)
    end
    m_plan = plan_rfft(m_pad, 1:3)
    field = create_zeros(3 * sim.n_total)
    energy = create_zeros(sim.n_total)
    return DemagPBC3D(nx, ny, nz, m_pad, M_pad, h_pad, tscale, m_plan, h_plan,
                      field, energy, "Demag")
end

#project the spectrum onto the longitudinal direction (positive kernel; the
#-1/mu0 in collect_h_energy turns it into the demag field):
#H_ij(k) = +k_i k_j/|k|^2 * prod(sinc^2(pi*f/n)) for k != 0, 0 for k = 0.
#kxu/kyu/kzu are 2*pi/(n*delta) so that k_i = f_i * k_iu; the cell shape factor
#only depends on the frequency indices (k_i*delta_i/2 = pi*f_i/n_i).
@kernel function spectral_project_kernel!(M_pad, nx::Int, ny::Int, nz::Int,
                                          kxu::T, kyu::T, kzu::T, invN::T) where {T<:AbstractFloat}
    a, b, c = @index(Global, NTuple)
    fx = T(a - 1)                                   #rfft half: non-negative freqs
    fy = (b - 1 <= ny ÷ 2) ? T(b - 1) : T(b - 1 - ny)
    fz = (c - 1 <= nz ÷ 2) ? T(c - 1) : T(c - 1 - nz)
    kx = fx * kxu
    ky = fy * kyu
    kz = fz * kzu
    k2 = kx * kx + ky * ky + kz * kz

    @inbounds Mx = M_pad[a, b, c, 1]
    @inbounds My = M_pad[a, b, c, 2]
    @inbounds Mz = M_pad[a, b, c, 3]

    if k2 == 0 #tin-foil convention: the uniform mode produces no field
        @inbounds M_pad[a, b, c, 1] = zero(T)
        @inbounds M_pad[a, b, c, 2] = zero(T)
        @inbounds M_pad[a, b, c, 3] = zero(T)
    else
        sx = fx == 0 ? one(T) : (sin(T(π) * fx / nx) / (T(π) * fx / nx))^2
        sy = fy == 0 ? one(T) : (sin(T(π) * fy / ny) / (T(π) * fy / ny))^2
        sz = fz == 0 ? one(T) : (sin(T(π) * fz / nz) / (T(π) * fz / nz))^2
        S = sx * sy * sz * invN
        D = (kx * Mx + ky * My + kz * Mz) / k2
        @inbounds M_pad[a, b, c, 1] = kx * D * S
        @inbounds M_pad[a, b, c, 2] = ky * D * S
        @inbounds M_pad[a, b, c, 3] = kz * D * S
    end
end

function effective_field(demag::DemagPBC3D, sim::MicroSim, spin::AbstractArray{T,1},
                         t::Float64; output=nothing) where {T<:AbstractFloat}
    mesh = sim.mesh
    nx, ny, nz = demag.nx, demag.ny, demag.nz

    distribute_m(spin, demag.m_pad, sim.mu0_Ms, nx, ny, nz)

    mul!(demag.M_pad, demag.m_plan, demag.m_pad)

    kernel! = spectral_project_kernel!(get_backend(spin), groupsize[])
    kernel!(demag.M_pad, nx, ny, nz,
            T(2π / (nx * mesh.dx)), T(2π / (ny * mesh.dy)), T(2π / (nz * mesh.dz)),
            demag.tscale; ndrange=(nx ÷ 2 + 1, ny, nz))

    inv_transform!(demag.h_pad, demag.M_pad, demag.h_plan, nx)

    heff = output == nothing ? demag.field : output

    collect_h_energy(heff, demag.energy, spin, demag.h_pad, sim.mu0_Ms,
                     T(mesh.volume), nx, ny, nz)

    return nothing
end
