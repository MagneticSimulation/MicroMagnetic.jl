using LinearAlgebra

# compute the demag field calculation using direct method
# the demag field is calculated as field = - demag_tensor*(spin*Ms)
# where demag_tensor is a NxN matrix with N the length of spin vector
mutable struct DirectDemag{T<:AbstractFloat} <: MicroEnergy
    tensor::AbstractArray{T,2}
    field::AbstractArray{T,1}
    energy::AbstractArray{T,1}
    name::String
end

function init_direct_demag(sim::MicroSim, Nx::Int, Ny::Int, Nz::Int)
    mesh = sim.mesh
    max_size = max(mesh.dx, mesh.dy, mesh.dz)
    dx = Float64(mesh.dx / max_size)
    dy = Float64(mesh.dy / max_size)
    dz = Float64(mesh.dz / max_size)

    nx = mesh.nx
    ny = mesh.ny
    nz = mesh.nz

    T = Float[]
    N = length(sim.spin)
    tensor = zeros(T, N, N)
    tensor_xx = zeros(T, nx, ny, nz)
    tensor_yy = zeros(T, nx, ny, nz)
    tensor_zz = zeros(T, nx, ny, nz)
    tensor_xy = zeros(T, nx, ny, nz)
    tensor_xz = zeros(T, nx, ny, nz)
    tensor_yz = zeros(T, nx, ny, nz)

    compute_demag_tensors(tensor_xx, tensors_kernel_xx!, Nx, Ny, Nz, dx, dy, dz)
    compute_demag_tensors(tensor_yy, tensors_kernel_yy!, Nx, Ny, Nz, dx, dy, dz)
    compute_demag_tensors(tensor_zz, tensors_kernel_zz!, Nx, Ny, Nz, dx, dy, dz)
    compute_demag_tensors(tensor_xy, tensors_kernel_xy!, Nx, Ny, Nz, dx, dy, dz)
    compute_demag_tensors(tensor_xz, tensors_kernel_xz!, Nx, Ny, Nz, dx, dy, dz)
    compute_demag_tensors(tensor_yz, tensors_kernel_yz!, Nx, Ny, Nz, dx, dy, dz)

    mu0_Ms = Array(sim.mu0_Ms)

    nxy = nx * ny
    for i in 1:nx, j in 1:ny, k in 1:nz
        id = nxy * (k - 1) + nx * (j - 1) + (i - 1)

        for ip in 1:nx, jp in 1:ny, kp in 1:nz
            idp = nxy * (kp - 1) + nx * (jp - 1) + (ip - 1)

            factor = -mu0_Ms[idp + 1] / mu_0

            di = abs(ip - i) + 1
            dj = abs(jp - j) + 1
            dk = abs(kp - k) + 1

            sx = ip - i >= 0 ? 1 : -1
            sy = jp - j >= 0 ? 1 : -1
            sz = kp - k >= 0 ? 1 : -1

            tensor[3 * id + 1, 3 * idp + 1] = factor * tensor_xx[di, dj, dk]       #xx
            tensor[3 * id + 1, 3 * idp + 2] = factor * sx * sy * tensor_xy[di, dj, dk] #xy
            tensor[3 * id + 1, 3 * idp + 3] = factor * sx * sz * tensor_xz[di, dj, dk] #xz

            tensor[3 * id + 2, 3 * idp + 1] = factor * sx * sy * tensor_xy[di, dj, dk] #yx = xy
            tensor[3 * id + 2, 3 * idp + 2] = factor * tensor_yy[di, dj, dk]       #yy
            tensor[3 * id + 2, 3 * idp + 3] = factor * sy * sz * tensor_yz[di, dj, dk] #yz

            tensor[3 * id + 3, 3 * idp + 1] = factor * sx * sz * tensor_xz[di, dj, dk]    #zx = xz
            tensor[3 * id + 3, 3 * idp + 2] = factor * sy * sz * tensor_yz[di, dj, dk]    #zy = yz
            tensor[3 * id + 3, 3 * idp + 3] = factor * tensor_zz[di, dj, dk]      #zz
        end
    end

    field = create_zeros(T, 3 * sim.n_total)
    energy = create_zeros(T, sim.n_total)

    demag = DirectDemag(kernel_array(tensor), field, energy, "Demag")
    return demag
end

function effective_field(demag::DirectDemag, sim::MicroSim, spin::AbstractArray{T,1},
                         t::Float64; output=nothing) where {T<:AbstractFloat}
    
    heff = output == nothing ? demag.field : output
    try
        mul!(heff, demag.tensor, spin)
    catch e #only works in cpu for AbstractFloat
        fill!(heff, 0)
        N = length(heff)
        for i = 1:N, j=1:N
            heff[i] += demag.tensor[i, j] * spin[j]
        end
    end

    N = sim.n_total
    factor = 0.5 * sim.mesh.volume
    #we borrow the zeeman kernel to compute the demag energy
    zeeman_kernel!(default_backend[])(spin, heff, demag.energy, sim.mu0_Ms, T(factor); ndrange=N)
    
    return nothing
end
